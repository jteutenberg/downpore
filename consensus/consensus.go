package consensus

import (
	"fmt"
	"github.com/jteutenberg/downpore/model"
	"github.com/jteutenberg/downpore/overlap"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/sequence/alignment"
)

//The consensus package contains a few utility functions for aligning reads based on matching seed sequences (overlaps) and
//generating consensus sequences from them
//The hard work is done in the alignment and overlap packages

func BuildConsensus(contig *overlap.SeedContig, sequences []sequence.Sequence, model model.Model, fullMatch bool) (*overlap.SeedContig, sequence.Sequence) {
	//fmt.Print("consensus from ",len(contig.Parts),"sequences...")
	k := 5
	costThreshold := uint(200) //get from model
	initialGapCost := uint(5)
	if model != nil {
		k = int(model.GetK())
		initialGapCost = uint(2)
	}
	seqs := make([][]uint16, 0, len(contig.Parts))
	rcs := make([]bool, 0, len(contig.Parts))
	seqMap := make([]int,0, len(contig.Parts)) //map used sequences to their indices in contig.Parts
	//pick contigs to use
	baseSeqIndex := -1
	for i, id := range contig.Parts {
		if contig.Matches != nil && contig.Matches[i].SeqA.GetID() == contig.Matches[i].SeqB.GetID() {
			baseSeqIndex = i
		}
		if contig.Approximate[i] {
			continue
		}
		b := sequences[id]
		start := contig.Offsets[i]
		if contig.Matches == nil {
			fmt.Println("Part",i," id/rc:",id,contig.ReverseComplement[i]," from ",contig.Offsets[i], " for ",contig.Lengths[i],"approx=",contig.Approximate[i])
		}
		if start < 0 {
			if start < -5 { //bad start, ignore
				continue
			}
			start = 0
		}
		end := contig.Offsets[i] + contig.Lengths[i]
		if end > b.Len() {
			if end > b.Len()+100 || (contig.ReverseComplement[i] && end > b.Len()+5) {
				continue //Ignore: bad end
			}
			end = b.Len()
		}
		if start >= end {
			fmt.Println("Nasty! The contig offset is past the end of the sequence!")
			start = end-1
		}
		b = b.SubSequence(start,end)
		if contig.ReverseComplement[i] {
			b = b.ReverseComplement()
		}

		rcs = append(rcs, contig.ReverseComplement[i])
		seqs = append(seqs, b.ShortKmers(k, false))
		seqMap = append(seqMap, i)
	}
	if len(seqs) < 3 {
		//fmt.Println("Down to",len(seqs),"from",len(contig.Parts),"due to approximate end matching.")
		return nil, nil
	}
	//fmt.Println("using",len(seqs),"ones with good start locations")

	//Perform consensus on these subsequences
	maxWarp := 16 //fixed to 2x16 actually
	var mod alignment.Measure
	if model != nil {
		mod = model.Clone()
	} else {
		mod = alignment.NewFivemerMeasure()
	}
	mod.SetSequences(seqs, rcs)

	dtw := alignment.NewDTWAligner(maxWarp,initialGapCost, mod, fullMatch, costThreshold, k)

	kmers,costs,pos := dtw.GlobalAlignment()
	ks := make([]uint16, 0, len(seqs[0]))
	var startPositions []int
	var endPositions []int
	for kmer := range kmers {
		<-costs
		endPositions = <-pos
		if startPositions == nil {
			startPositions = endPositions
		}
		ks = append(ks, kmer)
	}

	if len(ks) < 100 { //too short. Bad sequence match.
		fmt.Println("A bit short: ",len(ks))
		return nil, nil
	}

	consensusLen := len(ks)-k+1

	//update start and end positions of each sequence used
	for i := 0; i < len(contig.Lengths); i++ {
		contig.Lengths[i] = consensusLen
		contig.Approximate[i] = true
	}
	for i, index := range seqMap {
		contig.Approximate[index] = false
		if contig.ReverseComplement[index] {
			// in this case these values come from a reverse-complement of the sequence
			contig.Offsets[index] += len(seqs[i]) - endPositions[i] //offset is distance from the end
		} else {
			contig.Offsets[index] += startPositions[i]
		}
		contig.Lengths[index] = endPositions[i] - startPositions[i] + k - 1
	}
	var consensus sequence.Sequence
	if baseSeqIndex == -1 {
		consensus = sequence.NewByteSequenceFromKmers(-1,ks,k)
		//consensus = sequence.NewByteSubSequenceFromKmers(contig.Parts[ baseSeqIndex ],ks,k,consensusOffset,consensusInset)
	} else {
		//treat this as a subsequence, even though it is now has consensus contents
		consensusOffset := contig.Offsets[baseSeqIndex]
		consensusInset := contig.SeqLengths[baseSeqIndex] - consensusOffset - consensusLen
		consensus = sequence.NewByteSubSequenceFromKmers(contig.Parts[ baseSeqIndex ],ks,k,consensusOffset,consensusInset)
	}
	return contig, consensus
}
