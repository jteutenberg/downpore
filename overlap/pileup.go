package overlap

import (
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/util"
	"log"
	"math"
	"sort"
)

//functions in this file relate to multiple ordered overlaps between sequences that can be used to arrange the original sequences into a pileup
//It also provides functionality for pulling out multiple aligned sequences across the pileup


type Pileup struct {
	members []uint //sequences, ordered by start position
	starts []int
	ends []int
	//then useful reference points that are confident of sequence's relative positions
	referenceMembers [][]uint //index to members for each reference point
	referencePositions [][]uint //base position for each point given by referenceMembers
}

func (p *Pileup) Len() int {
	return len(p.members)
}
func (p *Pileup) Less(i, j int) bool {
	return p.starts[i] < p.starts[j]
}
func (p *Pileup) Swap(i, j int) {
	p.starts[i], p.starts[j] = p.starts[j], p.starts[i]
	p.ends[i], p.ends[j] = p.ends[j], p.ends[i]
	p.members[i], p.members[j] = p.members[j], p.members[i]
}

func (p *Pileup) MembersAt(offset int) []uint {
	ms := make([]uint, 0, 40)
	for i, start := range p.starts {
		if start >= offset {
			return ms
		}
		if p.ends[i] > offset {
			ms = append(ms,p.members[i])
		}
	}
	return ms
}

func (p *Pileup) MembersSpanning(from,to int) []uint {
	ms := make([]uint, 0, 40)
	for i, start := range p.starts {
		if start >= from {
			return ms
		}
		if p.ends[i] > to {
			ms = append(ms,p.members[i])
		}
	}
	return ms
}

func NewPileup(contigs []*SeedContig) *Pileup {
	//1. Get a list of all sequences

	maxID := 0
	for _, contig := range contigs {
		if contig == nil {
			continue
		}
		for _, p := range contig.Parts {
			if p > maxID {
				maxID = p
			}
		}
	}
	allSeqs := util.NewIntSetCapacity(maxID+1)
	for _, contig := range contigs {
		if contig == nil {
			continue
		}
		for _, p := range contig.Parts {
			allSeqs.Add(uint(p))
		}
	}
	pile := Pileup{members: allSeqs.AsUints()}

	backMap := make(map[uint]int)
	for i,m := range pile.members {
		backMap[m] = i
	}

	//2. Find a contig+offset start and end position for each sequence
	firstContig := make([]int, len(pile.members))
	lastContig := make([]int, len(pile.members))
	pile.starts = make([]int, len(pile.members)) //relative starts
	pile.ends = make([]int, len(pile.members)) //relative ends

	//estimated distance per contig
	contigOffsets := make([]int, len(contigs))
	seqEnds := make([]int, len(pile.members)) //final base of this sequence that is within a contig (or earliest base for rc sequences)
	for i,contig := range contigs {
		if contig == nil {
			if i > 0 {
				contigOffsets[i] = contigOffsets[i-1]+1000 //do something with these later.. should we remove its sequences??
			}
			continue
		}
		posEstimate := 0
		count := 0
		for j, p := range contig.Parts {
			rc := contig.ReverseComplement[j]
			index := backMap[uint(p)]
			if firstContig[index] == 0 {
				firstContig[index] = i
				//store the offset in starts for now
				if rc {
					pile.starts[index] = -( contig.SeqLengths[j] - (contig.Offsets[j]+contig.Lengths[j])) //we are *descending* through this sequence as we see more contigs
				} else {
					pile.starts[index] = -contig.Offsets[j] //the sequences typically overhang the front of the contig
				}
				if i == 0 && -pile.starts[index] > contigOffsets[0] {
					contigOffsets[0] = -pile.starts[index] //this way the earliest sequence will be given "start position" 0
				}
			}
			if i > 0 && lastContig[index] != 0 {
				//estimate the distance from the start if this contig to an earlier contig: its the difference in offsets
				//old := posEstimate
				if rc {
					posEstimate += contigOffsets[ lastContig[index] ] + contigs[lastContig[index]].Combined.Len() + seqEnds[index] - (contig.Offsets[j] + contig.Lengths[j])
					//log.Println(posEstimate/(count+1)," added ",posEstimate-old," m:",pile.members[index],"rc:",contig.ReverseComplement[j]," added from ",contig.Offsets[j]+contig.Lengths[j],"to",seqEnds[index]," that was ",i-lastContig[index],"contigs ago and ended at estimated",contigOffsets[ lastContig[index] ] + contigs[lastContig[index]].Combined.Len())
				} else {
					posEstimate += contigOffsets[lastContig[index] ] + contigs[lastContig[index]].Combined.Len() + contig.Offsets[j] - seqEnds[index]
					//log.Println(posEstimate/(count+1)," added",posEstimate-old," m:",pile.members[index],"rc:",contig.ReverseComplement[j]," added from ",seqEnds[index]," to ",contig.Offsets[j],"that was ",i-lastContig[index],"contigs ago (",lastContig[index],") and ended at estimated",contigOffsets[lastContig[index]],"+",contigs[lastContig[index]].Combined.Len()," = ",contigOffsets[ lastContig[index] ] + contigs[lastContig[index]].Combined.Len())
				}
				count++
			}
			lastContig[index] = i
			//and the relative base position of the end
			if rc {
				pile.ends[index] = contig.Combined.Len() + contig.Offsets[j] // easy for rc, as the offset gives the overhang
				seqEnds[index] = contig.Offsets[j] //position in the sequence (just the offset for rc sequences, as we count down)
			} else {
				//the "good" estimate of the combined length from the beginning of this contig , plus all the overhang at the end
				pile.ends[index] = contig.Combined.Len() + (contig.SeqLengths[j]-contig.Lengths[j]-contig.Offsets[j])
				seqEnds[index] = contig.Offsets[j] + contig.Lengths[j]
			}

		}
		if count > 0 {
			contigOffsets[i] = posEstimate/count
		} else if i > 0 {
			//oddly... a gap?! This contig must be garbage
			log.Println("Unable to estimate offset at ",i)
			contigOffsets[i] = contigOffsets[i-1]+1000 //do something with these later.. should we remove its sequences??
		}
	}
	//3. Assign final relative start/end points to each sequence, sort the sequences by start
	log.Println("Pileup of ",len(pile.members),"member sequences.")
	for index, _ := range pile.members {
		pile.starts[index] += contigOffsets[ firstContig[index] ]
		pile.ends[index] += contigOffsets[ lastContig[index] ]
	}
	sort.Sort(&pile) //by start

	log.Println("First seq is ",pile.members[0],"at",pile.starts[0],"to",pile.ends[0])
	log.Println("Last seq is ",pile.members[len(pile.members)-1],"at",pile.starts[len(pile.members)-1],"to",pile.ends[len(pile.members)-1])
	//we'll also make a reference point from most contigs
	pile.referenceMembers = make([][]uint,0,len(contigs))
	pile.referencePositions = make([][]uint,0,len(contigs))
	for i,contig := range contigs {
		if contig == nil {
			continue
		}
		//find all members at this point. Match new ones to the seed-space consensus?
		ms := pile.MembersAt(contigOffsets[i])
		log.Println("At contig ",i)
		log.Println("Have:",len(ms),"=",ms)
		log.Println("Expected:",len(contig.Parts),contig.Parts)
		//add a reference for the contig start

	}
	return &pile
}

//checkContainedSequence finds a subset of hits for this sequence that are roughly in the correct place relative to one another
func checkContainedSequence(id uint, futureContigs [][]*seeds.SeedMatch, seqSets []*util.IntSet, overlapSize, k int) (int, int){
	rightMost := len(futureContigs)-1
	for rightMost >= 1 && !seqSets[rightMost].Contains(id) {
		rightMost--
	}
	if rightMost == 0 {
		//only one hit. Just return.
		return 0,0
	}
	diagonal := make([]int, 0, rightMost+1)
	indices := make([]int, 0, rightMost+1)
	for i, set := range seqSets[:rightMost+1] {
		if set.Contains(id) {
			indices = append(indices, i)
			//find where this sequence is
			j := 0
			for j < len(futureContigs[i]) && futureContigs[i][j].SeqB.GetID() != int(id) {
				j++
			}
			match := futureContigs[i][j]
			//and determine the position on the diagonal
			if match.ReverseComplementQuery {
				diagonal = append(diagonal, match.SeqA.GetOffset()+match.SeqA.GetSeedOffset(match.MatchA[0],k) + match.SeqB.GetOffset() + match.SeqB.GetSeedOffset(match.MatchB[0],k))
			} else {
				diagonal = append(diagonal, match.SeqA.GetOffset()+match.SeqA.GetSeedOffset(match.MatchA[0],k) - match.SeqB.GetOffset() - match.SeqB.GetSeedOffset(match.MatchB[0],k))
			}
		}
	}
	//sort
	util.SortByValue(indices, diagonal)
	//and run a window across the diagonal, finding a good split to keep as many hits as possible
	window := overlapSize/2
	bestLength := 1
	bestStart := -1
	bestEnd := 0
	start := -1
	end := 0
	for start < len(indices)-bestLength {
		//step forward to the next spot
		start++
		first := diagonal[start]
		//move the end forward up to window length
		for end < len(indices) && first + window > diagonal[end] {
			end++
		}
		if end-start >= bestLength {
			bestLength = end-start
			bestStart = start
			bestEnd = end
		}
	}
	//remove all the others
	if bestLength == len(indices) {
		//all are good!
		return 0,rightMost
	} else if bestLength == 1 {
		//if we're dropping down to a single hit, just remove them all. Certain rubbish.
		bestLength = 0
	} else {
		for i := bestStart; i < bestEnd; i++ {
			diagonal[i] = indices[i]-math.MaxInt32
		}
		util.SortByValue(indices, diagonal)
	}
	for _, index := range indices[bestLength:] {
		set := seqSets[index]
		if set.Contains(id) {
			j := 0
			for j < len(futureContigs[index]) && futureContigs[index][j].SeqB.GetID() != int(id) {
				j++
			}
			if j < len(futureContigs[index])-1 {
				copy(futureContigs[index][j:], futureContigs[index][j+1:])
			}
			//futureContigs[index][j] = futureContigs[index][len(futureContigs[index])-1]
			futureContigs[index] = futureContigs[index][:len(futureContigs[index])-1]
			set.Remove(id)
		}
	}
	if bestLength == 0 {
		return -1,-1
	}
	return indices[0],indices[bestLength-1]
}

//hasOverhang checks whether the pileup (overlaps) match across a whole sequence or are overhanging, unmatched on one or both sides
func hasOverhang(id uint, overlaps [][]*seeds.SeedMatch, leftIndex, rightIndex, overlapSize, k int) (bool,int,int) {
	left := 0
	for left < len(overlaps[leftIndex]) && overlaps[leftIndex][left].SeqB.GetID() != int(id) {
		left++
	}
	//left is index in the first set containing id
	//right will be the index in the last set containing id
	var right int
	if leftIndex == rightIndex {
		right = left
	} else {
		for right < len(overlaps[rightIndex]) && overlaps[rightIndex][right].SeqB.GetID() != int(id) {
			right++
		}
	}
	leftMatch := overlaps[leftIndex][left]
	rightMatch := overlaps[rightIndex][right]
	var leftOverhang int
	var rightOverhang int
	if leftMatch.ReverseComplementQuery {
		//just switch left / right hits. They're in the correct direction individually.
		leftOverhang = rightMatch.SeqB.GetSeedOffset(rightMatch.MatchB[0],k)
		rightOverhang = leftMatch.SeqB.GetSeedOffsetFromEnd(leftMatch.MatchB[len(leftMatch.MatchB)-1],k)
	} else {
		leftOverhang = leftMatch.SeqB.GetSeedOffset(leftMatch.MatchB[0],k)
		rightOverhang = rightMatch.SeqB.GetSeedOffsetFromEnd(rightMatch.MatchB[len(rightMatch.MatchB)-1],k)
	}
	//start,_ := overlaps[leftIndex][left].GetAIndices(10)
	//_,end := overlaps[rightIndex][right].GetAIndices(10)
	//bSL , bEL := overlaps[leftIndex][left].GetBIndices(10)
	//bSR , bER := overlaps[rightIndex][right].GetBIndices(10)
	//fmt.Println("At",start,"to",end,"in big, and ",bSL,"to",bEL,"then",bSR,"to",bER,"in small RC=",leftMatch.ReverseComplementQuery,"/",rightMatch.ReverseComplementQuery,"Seq ",id," LHS overhang of ",leftOverhang, "from contig",leftIndex,"/",len(overlaps)," RHS overhang of ",rightOverhang," from contig ",rightIndex,"/",len(overlaps)," sequence length is",overlaps[leftIndex][left].SeqB.GetLength(),",",overlaps[rightIndex][right].SeqB.GetLength()," ids ",overlaps[leftIndex][left].SeqB.GetID(),"/",overlaps[rightIndex][right].SeqB.GetID(),"against limit of ",overlapSize*2)
	return (rightIndex < len(overlaps)-2 && rightOverhang > overlapSize*2) || (leftIndex > 1 && leftOverhang > overlapSize*2),left,right
}

//CleanupOverlaps removes all matches from the pileup (overlaps) which are not consistent, i.e. in-order and not overhanging too much
//It also returns the overlaps as a pileup
func CleanupOverlaps(overlaps [][]*seeds.SeedMatch, overlapSize, k int) {
	//1. make sets of each overlaps' sequences
	globalMax := 0
	seqSets := make([]*util.IntSet, len(overlaps))
	for i, overlap := range overlaps {
		max := 0
		for _, s := range overlap {
			id := s.SeqB.GetID()
			if id > max {
				max = id
			}
		}
		seqSets[i] = util.NewIntSetCapacity(max)
		for _, s := range overlap {
			seqSets[i].Add(uint(s.SeqB.GetID()))
		}
		if max > globalMax {
			globalMax = max
		}
	}

	checked := util.NewIntSetCapacity(globalMax+1)

	//2. left-to-right, find which sequences map as expected across the contigs. Discard those that don't belong.
	for i := 0; i < len(seqSets); i++ {
		seqs := seqSets[i]
		for ok,id := seqs.GetFirstID(); ok && !checked.Contains(id); ok,id = seqs.GetNextID(id) {
			//find the best contiguous match, and remove all others
			leftIndex, rightIndex := checkContainedSequence(id, overlaps[i:], seqSets[i:], overlapSize,k)
			if leftIndex == -1 {
				//in this case we removed all hits.
				continue
			}
			checked.Add(id) //either totally removed, or we are happy with its full span
			leftIndex += i
			rightIndex += i
			if overhangs, _,_ := hasOverhang(id, overlaps, leftIndex, rightIndex, overlapSize,k); overhangs {
				//TODO: add a flag. We can ignore those with overhangs if they match a reasonable chunk of the sequence if we like.
				for n := leftIndex; n <= rightIndex; n++ {
					if seqSets[n].Contains(id) {
						j := 0
						for j < len(overlaps[n]) && overlaps[n][j].SeqB.GetID() != int(id) {
							j++
						}
						if j < len(overlaps[n])-1 {
							copy(overlaps[n][j:],overlaps[n][j+1:])
						}
						//overlaps[n][j] = overlaps[n][len(overlaps[n])-1]
						overlaps[n] = overlaps[n][:len(overlaps[n])-1]
						seqSets[n].Remove(id)
					}
				}
			}
		}
	}
}
