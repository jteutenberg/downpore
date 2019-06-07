package commands

import (
	"bufio"
	"fmt"
	"github.com/jteutenberg/downpore/sequence"
	"log"
	"os"
	"strings"
)

type subseqCommand struct {
	args  map[string]string
	alias map[string]string
	desc  map[string]string
}

func NewSubSeqCommand() Command {
	args, alias, desc := MakeArgs(
		[]string{"input", "num_workers", "himem"},
		[]string{"", "4", "false"},
		[]string{"Fasta/fastq input file", "Number of worker threads to use", "Whether to cache reads in memory"})
	return &subseqCommand{args: args, alias: alias, desc: desc}
}
func (com *subseqCommand) GetName() string {
	return "subseq"
}
func (com *subseqCommand) GetArgs() (map[string]string, map[string]string, map[string]string) {
	return com.args, com.alias, com.desc
}

func (com *subseqCommand) Run(args map[string]string) {
	cacheReads := ParseBool(args["himem"])
	numWorkers := ParseInt(args["num_workers"])
	seqSet := sequence.NewFastaSequenceSet(args["input"], 0, numWorkers, cacheReads, true)

	ids := make(map[string]int) //map sequence names to ids in the sequence set
	seqs := seqSet.GetSequences()
	for seq := range seqs {
		name := seq.GetName()
		ids[name] = seq.GetID()
		if split := strings.Index(name," "); split > 0 && split < len(name) {
			ids[name[:split]] = seq.GetID()
		}
	}
	scanner := bufio.NewScanner(os.Stdin)
	for scanner.Scan() {
		line := scanner.Text()
		tokens := strings.Split(line, " ")

		start := ParseInt(tokens[0])
		end := ParseInt(tokens[1])
		rc := ParseBool(tokens[2])
		var name string
		if len(tokens) > 3 {
			name = tokens[3]
		}
		var seq sequence.Sequence
		if name != "" {
			if id, ok := ids[name]; ok {
				seqs := seqSet.GetNSequencesFrom(id, 1)
				seq = <-seqs
			} else {
				fmt.Println(name, "not found in", args["input"])
				log.Println(name, "not found in", args["input"])
			}
		} else {
			seqs := seqSet.GetNSequencesFrom(0, 1)
			seq = <-seqs
		}
		if seq != nil {
			if !strings.HasPrefix(seq.GetName(),name) {
				log.Println("Invalid name found:",seq.GetName()," != ",name)
				fmt.Println("Invalid name:",seq.GetName()," != ",name,"\n")
				continue
			}
			fmt.Printf(">%s_%d\n", seq.GetName(), start)
			if start > end { //circular wrapping
				sub1 := seq.SubSequence(start, seq.Len())
				sub2 := seq.SubSequence(0, end)
				if rc {
					fmt.Println(sub2.ReverseComplement().String() + sub1.ReverseComplement().String())
				} else {
					fmt.Println(sub1.String() + sub2.String())
				}
			} else {
				if end > seq.Len() {
					log.Println("Bad request:",start,"-",end," in sequence length ",seq.Len(),seq.GetName())
					end = seq.Len()
				}
				subseq := seq.SubSequence(start, end)
				if rc {
					fmt.Println(subseq.ReverseComplement())
				} else {
					fmt.Println(subseq)
				}
			}
		} else {
			fmt.Println("No sequence found.")
		}
	}
}
