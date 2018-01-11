
package commands

import (
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/trim"
	"log"
	"os"
)
type trimCommand struct {
	args map[string]string
	alias map[string]string
	desc map[string]string
}

func NewTrimCommand() Command {
	args, alias, desc := MakeArgs(
		[]string{"input", "end_threshold","middle_threshold","absolute_threshold","discard_middle","check_reads","adapter_threshold","extra_end_trim","extra_middle_trim","tag_adapters","verbosity","front_adapters","back_adapters","num_workers"},
		[]string{"","75","85","25","false","10000","90","5","100","true","1","","","4"},
		[]string{"Fasta/fastq/gzip input file","% identity for matching adapters at read edges","% identity for matching adapters that split reads","Identity (in bases) above which all adapters match","Whether to keep halves of split reads", "Number of reads to use to determine which adapters are present", "% identity required at check_adapters stage","Number of bases to remove around adapters at read edges", "Number of bases to remove around read-splitting adapters","Whether to add adapter names to output sequence names","Level (0-2) of output to stderr","Fasta/fastq file containing front adapters","Fasta/fastq file containing back adapters","Number of threads to use"})
	trim := trimCommand{args:args, alias:alias, desc:desc}
	return &trim
}
func (com *trimCommand) GetName() string {
	return "trim"
}

func (com *trimCommand) GetArgs() (map[string]string, map[string]string, map[string]string) {
	return com.args, com.alias, com.desc
}

func (com *trimCommand) Run(args map[string]string) {
	trimmer := trim.LoadTrimmer(args["front_adapters"],args["back_adapters"],6)
	seqSet := sequence.NewFastaSequenceSet(args["input"],50)
	numWorkers := ParseInt(args["num_workers"])
	trimmer.DetermineAdapters(seqSet, ParseInt(args["check_reads"]),ParseInt(args["adapter_threshold"]), numWorkers)
	trimmer.SetTrimParams(ParseInt(args["end_threshold"]),ParseInt(args["middle_threshold"]),ParseInt(args["absolute_threshold"]), ParseInt(args["extra_end_trim"]),ParseInt(args["extra_middle_trim"]), !ParseBool(args["discard_middle"]), ParseBool(args["tag_adapters"]))
	trimmer.Trim(seqSet, numWorkers)
	trimmer.PrintStats(seqSet)
	//and write
	log.Println("Writing trimmed sequences...")
	seqSet.Write(os.Stdout,true)
}
