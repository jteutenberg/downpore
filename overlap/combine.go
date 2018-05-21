package overlap

import (
	"github.com/jteutenberg/downpore/seeds"
)

type SeedContig struct {
	Combined          *seeds.SeedSequence //seed-based consensus
	Parts             []int
	ReverseComplement []bool
	Offsets           []int              //start position of Combined within each Part
	Lengths           []int              //length in bases of Part that contributes to Combined
	Approximate       []bool             //precision of the offset/length values
	SeqLengths        []int              //length in bases of the entire Part sequence
	Matches           []*seeds.SeedMatch //optional alignments of parts' seeds to the consensus
}

//gets the best seed as an index in shared SeqA, and the closest index in each SeqB
func trimToBestSeed(upto int, ms []*seeds.SeedMatch, k int) (consensus *seeds.SeedSequence, parts []*seeds.SeedSequence, cantTrim []bool) {
	parts = make([]*seeds.SeedSequence, len(ms), len(ms))
	cantTrim = make([]bool, len(ms), len(ms))
	//1. find the best front and back seeds
	bestCount := 0
	bestIndex := upto
	backCount := 0
	length := ms[0].SeqA.GetNumSeeds()
	backIndex := length - upto - 1
	for i := 0; i < upto; i++ {
		count := 0
		bCount := 0
		for _, match := range ms {
			for _, index := range match.MatchA {
				if index == i {
					count++
				}
				if index >= i {
					break
				}
			}
			for j := len(match.MatchA) - 1; j > 0; j-- {
				index := match.MatchA[j]
				if index == length-1-i {
					bCount++
				}
				if index <= length-1-i {
					break
				}
			}
		}
		if count-i > bestCount {
			bestCount = count - i //lose a bit for later seeds
			bestIndex = i
		}
		if bCount-i > backCount {
			backCount = bCount - i
			backIndex = length - 1 - i
		}
	}
	//2. Find the closest matches and trim as we go
	for j, match := range ms {
		//fmt.Println("Original offset is ",match.SeqB.GetOffset(),"and inset",match.SeqB.GetInset())
		index, bases, frontDistance := match.GetBaseIndex(bestIndex, k)     //this is an index + bases after. We trim by index + bases before...
		bIndex, backBases, backDistance := match.GetBaseIndex(backIndex, k) //This one is fine for trimming
		//if the final matching seed is before the back trim (or after the front trim)... ignore this sequence?
		if frontDistance > 50 || backDistance > 50 {
			parts[j] = match.SeqB
			cantTrim[j] = true
			//inaccurate positions
			continue
		}
		cantTrim[j] = false
		if bases > 0 && index < match.SeqB.GetNumSeeds()-1 {
			//move the seed forward one and the bases are those remaining in front
			bases = match.SeqB.GetNextSeedOffset(index, k) - bases
			index++
		} else if bases < 0 {
			bases = -bases
		}
		parts[j], _ = match.SeqB.Trimmed(bases, index, backBases, bIndex, k)
	}
	consensus, _ = ms[0].SeqA.Trimmed(0, bestIndex, 0, ms[0].SeqA.GetNumSeeds(), k)
	return consensus, parts, cantTrim
}

func NewSeedContig(ms []*seeds.SeedMatch, k int) *SeedContig {
	consensus, parts, trimFailed := trimToBestSeed(ms[0].SeqA.GetNumSeeds()/4, ms, k)

	contig := SeedContig{consensus, make([]int, len(ms), len(ms)), make([]bool, len(ms), len(ms)), make([]int, len(ms), len(ms)), make([]int, len(ms), len(ms)), trimFailed, make([]int, len(ms), len(ms)), ms}
	for i, part := range parts {
		contig.Parts[i] = part.GetID()
		contig.ReverseComplement[i] = part.IsReverseComplement()
		parent := part
		for parent.Parent != nil {
			parent = parent.Parent
		}
		contig.SeqLengths[i] = parent.Len()
		contig.Offsets[i] = part.GetOffset()
		contig.Lengths[i] = parent.Len() - part.GetOffset() - part.GetInset()
	}
	return &contig
}

//Remove leaves the combined seed sequence but removes the given part and all its related information
func (contig *SeedContig) Remove(part int) {
	index := 0
	for ;index < len(contig.Parts); index++ {
		if contig.Parts[index] == part {
			break
		}
	}
	last := len(contig.Parts)-1
	if last != index {
		contig.Parts[index] = contig.Parts[last]
		contig.Lengths[index] = contig.Lengths[last]
		contig.Offsets[index] = contig.Offsets[last]
		contig.SeqLengths[index] = contig.SeqLengths[last]
		contig.ReverseComplement[index] = contig.ReverseComplement[last]
		contig.Approximate[index] = contig.Approximate[last]
		contig.Matches[index] = contig.Matches[last]
		contig.Matches[last] = nil
	}
	contig.Parts = contig.Parts[:last]
	contig.Lengths = contig.Lengths[:last]
	contig.Offsets = contig.Offsets[:last]
	contig.SeqLengths = contig.SeqLengths[:last]
	contig.ReverseComplement = contig.ReverseComplement[:last]
	contig.Approximate = contig.Approximate[:last]
	contig.Matches = contig.Matches[:last]
}

func BuildConsensus(sg *seeds.SeedIndex, overlaps []*seeds.SeedMatch) *SeedContig {
	k := int(sg.GetSeedLength())
	anchors := make([]int, 0, len(overlaps)+1) //start seed indices
	anchorOffsets := make([]int, 0, len(overlaps)+1)
	mismatches := make([]int, 0, len(overlaps)+1)
	seqs := make([]*seeds.SeedSequence, 0, len(overlaps)+1)

	for _, lap := range overlaps {
		if lap.ReverseComplementQuery {
			lap.ReverseComplement(k,sg)
			//This means SeqA is identical (i.e. by reference) to the rest
			//But SeqB is a new SeedSequence with reversed edge offsets, etc. <-- it no longer relates to the read
		}
	}
	//generate anchor points by estimating the location of the middle seed
	fullSeeds := overlaps[0].SeqA.GetNumSeeds()
	mid := fullSeeds / 2
	for _, lap := range overlaps {
		anchor, offset, _ := lap.GetBaseIndex(mid, k)
		s := lap.SeqB
		if anchor >= s.GetNumSeeds() || len(lap.MatchA) < 3 { //1-2 seeds.. just not enough
			continue
		}
		ca, cb := lap.GetBasesCovered(k)
		if ca < 25 || cb < 25 { //~30 bases minimum too. Can occur when seeds overlap
			continue
		}

		//and now trim to just the overlap
		s, trimmed := s.Trimmed(overlaps[0].SeqA.GetSeedOffset(lap.MatchA[0], k), lap.MatchB[0], overlaps[0].SeqA.GetSeedOffsetFromEnd(lap.MatchA[len(lap.MatchA)-1], k), lap.MatchB[len(lap.MatchB)-1], k)
		anchor -= trimmed
		if anchor < 0 || anchor > s.GetNumSeeds() {
			//in this case, s only partially overlaps the query, missing the crucial centre seed
			continue
		}
		anchors = append(anchors, anchor)
		anchorOffsets = append(anchorOffsets, offset)
		mismatches = append(mismatches, fullSeeds-len(lap.MatchA))
		seqs = append(seqs, s)
	}
	if len(seqs) > 1 {
		overlap := seeds.Consensus(seqs, mismatches, anchors, anchorOffsets, k)
		if len(overlap) > 1 {
			return NewSeedContig(overlap, k)
		}
	}
	return nil
}
