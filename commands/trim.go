package commands

import (
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/trim"
	"log"
	"os"
)

type trimCommand struct {
	args  map[string]string
	alias map[string]string
	desc  map[string]string
}

func NewTrimCommand() Command {
	args, alias, desc := MakeArgs(
		[]string{"input", "k", "chunk_size", "middle_threshold", "discard_middle", "check_reads", "adapter_threshold", "extra_end_trim", "extra_middle_trim", "tag_adapters", "verbosity", "front_adapters", "back_adapters", "num_workers", "himem"},
		[]string{"", "6", "5000", "85", "false", "10000", "90", "5", "100", "true", "1", "", "", "4", "false"},
		[]string{"Fasta/fastq/gzip input file", "k-mer size to use when matching adapters", "Split long reads into chunks of this size when indexing", "% identity for matching adapters that split reads", "Whether to keep halves of split reads", "Number of reads to use to determine which adapters are present", "% identity required at check_adapters stage", "Number of bases to remove around adapters at read edges", "Number of bases to remove around read-splitting adapters", "Whether to add adapter names to output sequence names", "Level (0-2) of output to stderr", "Fasta/fastq file containing front adapters", "Fasta/fastq file containing back adapters", "Number of threads to use", "Whether to cache all reads in memory"})
	trim := trimCommand{args: args, alias: alias, desc: desc}
	return &trim
}
func (com *trimCommand) GetName() string {
	return "trim"
}

func (com *trimCommand) GetArgs() (map[string]string, map[string]string, map[string]string) {
	return com.args, com.alias, com.desc
}

func (com *trimCommand) Run(args map[string]string) {
	numWorkers := ParseInt(args["num_workers"])
	trimmer := trim.LoadTrimmer(args["front_adapters"], args["back_adapters"], ParseInt(args["k"]))
	seqSet := sequence.NewFastaSequenceSet(args["input"], 50, numWorkers, ParseBool(args["himem"]), false)
	trimmer.SetVerbosity(ParseInt(args["verbosity"]))
	trimmer.DetermineAdapters(seqSet, ParseInt(args["check_reads"]), ParseInt(args["adapter_threshold"]), numWorkers)
	trimmer.SetTrimParams(ParseInt(args["middle_threshold"]), ParseInt(args["extra_end_trim"]), ParseInt(args["extra_middle_trim"]), ParseInt(args["chunk_size"]), !ParseBool(args["discard_middle"]), ParseBool(args["tag_adapters"]))
	trimmer.Trim(seqSet, numWorkers)
	trimmer.PrintStats(seqSet)
	//and write
	log.Println("Writing trimmed sequences...")
	seqSet.Write(os.Stdout, true)
}
