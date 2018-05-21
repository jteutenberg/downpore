package commands

import (
	"fmt"
	"github.com/jteutenberg/downpore/model"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/sequence/alignment"
)

type alignCommand struct {
	args  map[string]string
	alias map[string]string
	desc  map[string]string
}

func NewAlignCommand() Command {
	args, alias, desc := MakeArgs(
		[]string{"input", "rc_input", "model", "reference"},
		[]string{"", "", "", ""},
		[]string{"Fasta/fastq input file", "Additional input file containing sequences from reverse-complement reads", "Model file containing current levels", "(optional) A fasta file containing a reference sequence to align against"})
	cons := alignCommand{args: args, alias: alias, desc: desc}
	return &cons
}
func (com *alignCommand) GetName() string {
	return "align"
}

func (com *alignCommand) GetArgs() (map[string]string, map[string]string, map[string]string) {
	return com.args, com.alias, com.desc
}

func (com *alignCommand) Run(args map[string]string) {
	var m alignment.Measure
	k := 5
	initialGapCost := uint(5)
	costThreshold := uint(200)
	if len(args["model"]) > 0 {
		mod := model.NewModel(args["model"], false)
		k = int(mod.GetK())
		m = mod
		costThreshold = uint(200) //get from model
		initialGapCost = uint(2)
	} else {
		m = alignment.NewFivemerMeasure()
	}
	//read each sequence, convert to short kmers
	kmerSeqs := make([][]uint16, 0, 100)
	seqSet := sequence.NewFastaSequenceSet(args["input"], 0, 1, false, false)
	seqs := seqSet.GetSequences()
	for seq := range seqs {
		kmerSeqs = append(kmerSeqs, seq.ShortKmers(k, false))
	}
	nonRC := len(kmerSeqs)
	seqSet = sequence.NewFastaSequenceSet(args["rc_input"], 0, 1, false, false)
	seqs = seqSet.GetSequences()
	for seq := range seqs {
		kmerSeqs = append(kmerSeqs, seq.ShortKmers(k, false))
	}

	var ref []uint16
	if args["reference"] != "" {
		seqSet = sequence.NewFastaSequenceSet(args["reference"], 0,1, false, false)
		seqs = seqSet.GetSequences()
		seq := <-seqs
		ref = seq.ShortKmers(k,false)
		for _, ok := <-seqs; ok; _, ok = <-seqs {
		} //drain the channel
	}

	maxWarp := 16

	dtw := alignment.NewDTWAligner(maxWarp, initialGapCost, m, false , costThreshold, k)
	//align 'em
	rc := make([]bool, len(kmerSeqs), len(kmerSeqs))
	for ; nonRC < len(rc); nonRC++ {
		rc[nonRC] = true
	}
	m.SetSequences(kmerSeqs, rc)
	var kmers <-chan uint16
	var costs <-chan *alignment.QualityMetrics
	var positions <-chan []int
	if ref == nil {
		kmers, costs, positions = dtw.GlobalAlignment()
	} else {
		kmers, costs, positions = dtw.GlobalAlignmentTo(ref)
	}

	finalCosts := make([]uint, len(kmerSeqs), len(kmerSeqs))
	prevPos := make([]int, len(kmerSeqs), len(kmerSeqs))
	lines := make([]string, len(kmerSeqs)+1, len(kmerSeqs)+1)
	for kmer := range kmers {
		ks := sequence.KmerString(int(kmer), k)
		mid := ks[len(ks)/2 : len(ks)/2+1]
		cs := <-costs
		pos := <-positions
		skips := 1
		for i, p := range pos {
			sk := p - prevPos[i]
			if sk > skips {
				skips = sk
			}
			finalCosts[i] = cs.CostDelta
		}
		// all lines will get exactly 'skip' new characters added
		for i := 1; i < skips; i++ {
			lines[0] += "." // a gap
		}
		lines[0] += mid
		//consensus has been written, now add the sequences
		for i, p := range pos {
			//fmt.Println("seq",i,"position",p,"from",prevPos[i])
			sk := p - prevPos[i] //the skip for this sequence
			if sk <= 0 {         //at worst a stay. Otherwise we're a bit mixed up
				for j := 0; j < skips; j++ {
					lines[i+1] += "."
				}
				continue
			}
			bases := skips //how many to write out
			nextKmer := sequence.KmerString(int(kmerSeqs[i][p]), k)
			//write out skips, but save some to just write out as kmers (usually we write nothing here)
			for ; sk > len(nextKmer)/2+1; sk-- {
				oldMer := sequence.KmerString(int(kmerSeqs[i][p-sk]), k)[len(nextKmer)/2 : len(nextKmer)/2+1]
				lines[i+1] += oldMer
				//lines[i+1] += "." //for large skips (only initial gap?)
				bases--
			}
			//add at least one base, or more if there are skips involved
			mid := nextKmer[len(nextKmer)/2+1-sk : len(nextKmer)/2+1] //remaining sk as bases
			bases -= len(mid)
			for ; bases > 0; bases-- { //and fill in gaps against larger skips. Umm.. should be none?
				lines[i+1] += "."
			}
			//put in the bases we left aside
			lines[i+1] += mid
		}
		prevPos = pos
	}
	for _, line := range lines {
		fmt.Println(line)
	}
}
