package commands

import (
	"fmt"
	"github.com/jteutenberg/downpore/model"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/sequence/alignment"
	"github.com/jteutenberg/downpore/util/sequtil"
)

type consensusCommand struct {
	args  map[string]string
	alias map[string]string
	desc  map[string]string
}

func NewConsensusCommand() Command {
	args, alias, desc := MakeArgs(
		[]string{"input", "rc_input", "model","matrix","k"},
		[]string{"", "", "","","5"},
		[]string{"Fasta/fastq input file", "Additional input file containing sequences from reverse-complement reads", "Model file containing current levels","K-mer confusion matrix to use in place of a model","K-mer size for alignment when no model specified"})
	cons := consensusCommand{args: args, alias: alias, desc: desc}
	return &cons
}
func (com *consensusCommand) GetName() string {
	return "consensus"
}

func (com *consensusCommand) GetArgs() (map[string]string, map[string]string, map[string]string) {
	return com.args, com.alias, com.desc
}

func (com *consensusCommand) Run(args map[string]string) {
	var m alignment.Measure
	k := ParseInt(args["k"])
	maxWarp := 16
	initialGapCost := uint(5)
	costThreshold := uint(200)
	if len(args["model"]) > 0 {
		mod := model.NewModel(args["model"], false)
		k = int(mod.GetK())
		m = mod
		costThreshold = uint(200) //get from model
		initialGapCost = uint(2)
	} else if len(args["matrix"]) > 0 {
		matrix, mk := sequtil.LoadConfusionMatrix(args["matrix"])
		k = mk
		m = alignment.NewMatrixDistance(k,matrix)
	} else {
		if k == 1 {
			m = alignment.NewBaseMeasure()
		} else if k <= 3 {
			k = 3
			m = alignment.NewThreemerMeasure()
		} else if k == 4 {
			m = alignment.NewFourmerMeasure()
		} else if k == 5 {
			m = alignment.NewFivemerMeasure()
		} else {
			//k = 6
			//m = alignment.NewSixmerMeasure()
			k = 5
			m = alignment.NewEditDistance(k,3,4,1)
		}
	}

	//read each sequence, convert to short kmers
	kmerSeqs := make([][]uint16, 0, 100)
	seqSet := sequence.NewFastaSequenceSet(args["input"], 0, 1, false, false)
	seqs := seqSet.GetSequences()
	for seq := range seqs {
		kmerSeqs = append(kmerSeqs, seq.ShortKmers(k,false))
	}
	nonRC := len(kmerSeqs)
	seqSet = sequence.NewFastaSequenceSet(args["rc_input"], 0, 1, false, false)
	seqs = seqSet.GetSequences()
	for seq := range seqs {
		kmerSeqs = append(kmerSeqs, seq.ShortKmers(k,false))
	}
	dtw := alignment.NewDTWAligner(maxWarp, initialGapCost, m, false, costThreshold, k)
	//dtw := alignment.NewDTWAligner(maxWarp, initialGapCost, m, true, costThreshold, k)
	//consensus 'em
	rc := make([]bool, len(kmerSeqs))
	for ; nonRC < len(rc); nonRC++ {
		rc[nonRC] = true
	}
	m.SetSequences(kmerSeqs, rc)
	kmers, costs, finalResult := dtw.GlobalConsensus()
	first := true
	var costsString string //three quality metrics
	for i := k; i > 0; i-- {
		costsString += "."
	}
	votesString := costsString
	spaceString := costsString

	consensus := make([]uint16, 0, 1000)
	for kmer := range kmers {
		cost := <-costs
		dc := cost.CostDelta
		if dc > 0 {
			dc = 1 + dc/30
			if dc >= 10 {
				dc = 9
			}
		}
		sp := int(cost.StateSpaceSize / 2)
		if sp > 7 {
			if sp > 50 {
				sp = 9
			} else {
				sp = 8
			}
		}
		costsString = fmt.Sprint(costsString, dc)
		votesString = fmt.Sprint(votesString, int(cost.ExactFraction*9.99))
		spaceString = fmt.Sprint(spaceString, sp)
		consensus = append(consensus, kmer)
		if first {
			fmt.Print(sequence.KmerString(int(kmer), k))
			first = false
		} else {
			ks := sequence.KmerString(int(kmer), k)
			fmt.Print(string(ks[len(ks)-1]))
		}
	}
	fmt.Println()
	fmt.Println(costsString)
	fmt.Println(votesString)
	fmt.Println(spaceString)
	<-finalResult
	//fmt.Println("\n",cost,"over",len(consensus),"bases",dtw.ConsensusCost(consensus))
}
