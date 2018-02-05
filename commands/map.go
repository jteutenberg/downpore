package commands

import (
	"fmt"
	"github.com/jteutenberg/downpore/mapping"
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/util/sequtil"
	"log"
)

type mapCommand struct {
	args  map[string]string
	alias map[string]string
	desc  map[string]string
}

func NewMapCommand() Command {
	args, alias, desc := MakeArgs(
		[]string{"input", "reference"},
		[]string{"", ""},
		[]string{"Fasta/fastq input file", "A fasta file containing a reference sequence to align against"})
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
	fmt.Println("reference length", reference.Len())
	k := 10

	kmerCounts := sequtil.KmerOccurrences(seqSet.GetSequences(), k, 1)
	values := make([]float64, len(kmerCounts))
	for i, count := range kmerCounts {
		j := seeds.ReverseComplement(uint(i), uint(k))
		ratio := float64(count) / float64(kmerCounts[j]+1)
		if ratio < 1.0 {
			ratio = 1.0 / ratio
		}
		values[i] = ratio //replacing values
	}

	log.Println("Counting complete. Preparing to start indexing and querying...")
	bottom, top := sequtil.TopOccurrences(kmerCounts, uint(k), len(kmerCounts)/100, len(kmerCounts)/50)
	for _, x := range bottom {
		values[x] = 10000
	}
	for _, x := range top {
		values[x] = 10000
	}

	mapper := mapping.NewMapper(reference, uint(k), values, 75, 4)
	//read each sequence, map against reference
	seqSet = sequence.NewFastaSequenceSet(args["input"], 50, 1, false, false)
	seqs := seqSet.GetSequences()
	for seq := range seqs {
		maps := mapper.Map(seq)
		fmt.Print(seqSet.GetName(seq.GetID()), " ")
		for _, m := range maps {
			fmt.Print(m.Start, "-", m.End, "(", m.RC, ") ")
		}
		fmt.Println()
	}
}
