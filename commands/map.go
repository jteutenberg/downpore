package commands

import (
	"fmt"
	"github.com/jteutenberg/downpore/mapping"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/util/sequtil"
	"log"
	"os"
	"runtime/pprof"
)

type mapCommand struct {
	args  map[string]string
	alias map[string]string
	desc  map[string]string
}

func NewMapCommand() Command {
	args, alias, desc := MakeArgs(
		[]string{"input", "reference","circular","k","num_workers"},
		[]string{"", "","true","10","4"},
		[]string{"Fasta/fastq input file", "A fasta file containing a reference sequence to align against","Whether the reference genome is circular","Length of seeds in bases","The number of worker process to use for mapping"})
	cons := mapCommand{args: args, alias: alias, desc: desc}
	return &cons
}
func (com *mapCommand) GetName() string {
	return "map"
}

func (com *mapCommand) GetArgs() (map[string]string, map[string]string, map[string]string) {
	return com.args, com.alias, com.desc
}

func (com *mapCommand) Run(args map[string]string) {
	seqSet := sequence.NewFastaSequenceSet(args["reference"], 0, 1, false, false)
	refs := seqSet.GetSequences()
	reference := <-refs
	k := ParseInt(args["k"])
	numWorkers := ParseInt(args["num_workers"])
	minLength := 500
	circular := ParseBool(args["circular"])
	overlap := 1000

	f, _ := os.Create("./cprof")
	pprof.StartCPUProfile(f)
	kmerCounts := sequtil.KmerOccurrences(seqSet.GetSequences(), k, numWorkers)
	values := make([]float64, len(kmerCounts))
	var tot uint64
	for _, count := range kmerCounts {
		tot += count
	}
	tf := float64(tot)
	//aim for fairly low frequency: about 1:200000 bases?
	targetFreq := 0.000005
	for i, count := range kmerCounts {
		freq := float64(count)/tf
		if freq <= targetFreq {
			values[i] = targetFreq - freq
		} else {
			values[i] = freq - targetFreq
		}
	}

	log.Println("K-mer counting complete. Preparing to start indexing and querying...")
	bottom, top := sequtil.TopOccurrences(kmerCounts, uint(k), len(kmerCounts)/100, len(kmerCounts)/50)
	for _, x := range bottom {
		values[x] = 10000
	}
	for _, x := range top {
		values[x] = 10000
	}

	mapper := mapping.NewMapper(reference, circular, uint(k), values, 75, overlap, 4)
	//read each sequence, map against reference
	seqSet = sequence.NewFastaSequenceSet(args["input"], minLength, 1, false, false)
	seqs := seqSet.GetSequences()
	unmapped := 0
	mapped := 0
	multiple := 0
	total := 0

	results := make(chan []*mapping.Mapping, numWorkers*2)
	done := make(chan bool, numWorkers)
	for i := 0; i < numWorkers; i++ {
		go mapper.MapWorker(seqs, results, done)
	}
	allDone := make(chan bool, 1)
	go func() {
		for maps := range results {
			if len(maps) > 0 {
				/*fmt.Print(seqSet.GetName(maps[0].Query.GetID()), " ")
				for _, m := range maps {
					fmt.Print(m.Start, "-", m.End, "(", m.RC, ") ")
				}
				fmt.Println(" Len:",maps[0].Query.Len())*/
				for _, m := range maps {
					fmt.Println(mapper.AsString(m))
				}
				if len(maps) == 1 {
					mapped++
				} else {
					multiple++
				}
				total += len(maps)
			} else {
				unmapped++
			}
		}
		allDone <- true
	}()

	for i := 0; i < numWorkers; i++ {
		<-done
	}
	close(results)
	<-allDone
	os.Stderr.WriteString(fmt.Sprintln("Uniquely mapped:",mapped))
	os.Stderr.WriteString(fmt.Sprintln("Multiple mappings:",multiple))
	os.Stderr.WriteString(fmt.Sprintln("total:",total))
	os.Stderr.WriteString(fmt.Sprintln("Unmapped:",unmapped))
	pprof.StopCPUProfile()
	f.Close()
}

