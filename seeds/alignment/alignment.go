package alignment

import (
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/util"
	"log"
)

type Aligner interface {
	PairwiseAlignments(*seeds.SeedSequence, *seeds.SeedSequence, *util.IntSet, *util.IntSet, int,int,bool) []*seeds.SeedMatch
	//PairwiseBestAlignment([]int,[]int, *util.IntSet, *util.IntSet, int,int) *seeds.SeedMatch
	//GlobalConsensus(*[]seeds.SeedSequence) []*seeds.SeedMatch
}

type seedAligner struct {
	stackPool []pairState //for reuse
	statesStack []*pairState //pointers to the stack pool
	nextState int //index of top of statesStack

	reduced []int //temporary reduced sequence
	aMapping []int //temporary mapping from reduced to original indices
	open []*pairState //pointing to allStates
	initials []*pairState //states with no matches yet
	results []*pairState
}

type pairState struct {
	aPos int //position of last matching seeds
	bPos int
	aGap int //number of bases since last aPos
	bGap int //number of bases since last bPos
	aGapIndex int //position of seed at gap aGap from aPos
	length int //length of the prev chain in seeds
	prev *pairState
	stackIndex int
}

//NewSeedAligner maxLength is the longest input sequence (for first argument to an alignment function)
func NewSeedAligner(maxLength int) *seedAligner {
	sa := seedAligner{stackPool: make([]pairState, 10000), statesStack: make([]*pairState,10000), reduced: make([]int,maxLength), open:make([]*pairState, 500), initials:make([]*pairState, maxLength), results:make([]*pairState, 500)}
	sa.nextState = len(sa.statesStack)-1
	sa.aMapping = make([]int, len(sa.reduced)/2)
	for i, _ := range sa.stackPool {
		sa.statesStack[i] = &sa.stackPool[i]
	}
	return &sa
}

func (align *seedAligner) popState() *pairState {
	s := align.statesStack[align.nextState]
	s.stackIndex = align.nextState
	align.nextState--
	return s
}

func (align *seedAligner) pushState(s *pairState) {
	//find where this is and swap with one above the stack top
	n := align.nextState+1
	top := align.statesStack[n]
	align.statesStack[n] = s
	align.statesStack[s.stackIndex] = top
	top.stackIndex = s.stackIndex
	s.stackIndex = n
	align.nextState = n
}

func (s *pairState) extractMatch(aMapping []int) *seeds.SeedMatch {
	ma := make([]int,s.length)
	mb := make([]int,s.length)
	for s != nil {
		ma[s.length-1] = aMapping[s.aPos/2]
		mb[s.length-1] = s.bPos/2
		s = s.prev
	}
	return &seeds.SeedMatch{MatchA:ma, MatchB:mb}
}

func (align *seedAligner) reset() {
	align.nextState = len(align.statesStack)-1
}

func (align *seedAligner) prepareInitial(aSegments []int, bSet *util.IntSet, minMatches, k int) (int,[]int,[]int) {
	maxAIndex := len(aSegments)-minMatches*2+1
	//fill the initial states from a
	aLen := 0
	offset := -k
	startSize := 0
	aRed := align.reduced
	aMapping := align.aMapping
	prevSeed := -1
	for i := 1; i < len(aSegments); i += 2 {
		aSeed := aSegments[i]
		if !bSet.Contains(uint(aSeed)) {
			offset += aSegments[i-1] + k
			maxAIndex--
			continue
		}
		if aSeed == prevSeed && (i >= len(aSegments)-2 || aSegments[i+2] == prevSeed) {
			offset += aSegments[i-1] + k
			maxAIndex--
			continue
		}
		prevSeed = aSeed

		offset += aSegments[i-1] + k
		aRed[aLen*2] = offset
		aRed[aLen*2+1] = aSeed
		aMapping[aLen] = i/2
		offset = -k
		if aLen <= maxAIndex {
			//set up an initial state
			state := align.popState()
			state.aPos = aLen*2+1 //position in reduced sequence
			state.length = 0 //unmatched so far
			state.prev = nil
			align.initials[aLen] = state
			startSize++
		}
		aLen++
	}
	aRed[aLen*2] = 0

	//trim back from startSize for any states past the new maxAIndex
	for startSize > 0 && align.initials[startSize-1].aPos > maxAIndex {
		startSize--
		align.pushState(align.initials[startSize])
	}
	return startSize, aRed[:aLen*2+1], aMapping[:aLen]
}

func (align *seedAligner) removeOpenState(index, minMatches, openSize, resultsSize int) (int,int,int) {
	s := align.open[index]
	align.open[index] = align.open[openSize-1]
	openSize--
	if s.length >= minMatches {
		if (s.length*2)/3 > minMatches {
			minMatches = (s.length*2)/3
		}
		align.results[resultsSize] = s
		resultsSize++

	} else {
		//discard the entire chain
		for s != nil {
			align.pushState(s)
			s = s.prev
		}
	}
	return openSize, resultsSize, minMatches
}

func (align *seedAligner) PairwiseAlignments(a, b *seeds.SeedSequence, aSet, bSet *util.IntSet, minMatches, k int, debug bool) []*seeds.SeedMatch {
	aSegments := a.GetSegments()
	bSegments := b.GetSegments()
	//effectively a slice into allStates
	if minMatches == 0 {
		minMatches = 1
	}
	align.reset()
	initialSize, aRed, aMapping := align.prepareInitial(aSegments, bSet, minMatches, k)
	//log.Println(aRed)
	
	openSize := 0
	resultsSize := 0
	//TODO ? two-pass generate counts of each seed in the sequences
	// - decrement b's count as we traverse it. If seed reaches 0, we can remove those initial states.

	//step through b
	bLen := len(bSegments)
	maxBIndex := len(bSegments)-minMatches*2+1
	bOffset := 0 //extra offset (gap) to the next bSeed
	prevSeed := -1
	for bIndex := 1; bIndex < bLen; bIndex+=2 {
		bSeed := bSegments[bIndex]
		if !aSet.Contains(uint(bSeed)) {
			bOffset += bSegments[bIndex+1] + k
			continue
		}
		if bSeed == prevSeed && (bIndex >= len(bSegments)-2 || bSegments[bIndex+2] == prevSeed) {
			bOffset += bSegments[bIndex+1] + k
			continue
		}
		prevSeed = bSeed
		if debug {
			log.Println("=============",bSeed,"(",bIndex/2,")===========")
		}
		found := -1 //first aIndex matched
		prevFound := -1 //if found != -1, this is the "i" (aIndex) at which it was set. This is a *possible* matching index, as the open list does get shuffled occassionally.
		//for each open state, can we find a suitable match for bSeed in reduced aSegments?
		searchMatch:
		for i := openSize-1; i >= 0; i-- {
			s := align.open[i]
			s.bGap += bOffset
			//look for any good match close after the current aPos
			minGap := (s.bGap*2)/3 - k
			maxGap := (s.bGap*3)/2 + k +1
			if minGap < 0 {
				minGap = -k
				if maxGap < 0 {
					maxGap = 0
				}
			} else if maxGap < 20 {
				maxGap = 20
				minGap = 0
			}
			if debug {
				log.Println("checking gap between",minGap,maxGap,"from index",s.aPos,"(gap index",s.aGapIndex,")")
			}
			//walk forward the aGap until it is within range
			ended := false
			for s.aGap < minGap {
				//if this has hit the end, add it to results (if long enough)
				if s.aGapIndex >= len(aRed) {
					//log.Println("state ended when searching at",s.aGapIndex)
					ended = true
					openSize, resultsSize, minMatches = align.removeOpenState(i,minMatches, openSize, resultsSize)
					break searchMatch
				}
				s.aGap += aRed[s.aGapIndex+1] + k
				s.aGapIndex += 2
			}
			//check all seeds in this narrow band
			if !ended {
				if debug {
					log.Println("Ended up with a gap at",s.aGap,"having walked to index",s.aGapIndex)
				}
				if s.aGap <= maxGap {
					g := s.aGap
					for j := s.aGapIndex; j < len(aRed) && g <= maxGap; j+=2 {
						if aRed[j] == bSeed {
							if debug {
								log.Println("matched.")
							}
							//if found is already set ==> found 2+. Check for redundancy / dominated chains
							if found != -1 && prevFound > i && prevFound < openSize {
								s2 := align.open[prevFound]
								if s.aPos == s2.aPos && s.bPos == s2.bPos {
									//which to remove?
									if s.length < s2.length {
										openSize, _, _ = align.removeOpenState(i, openSize, resultsSize, s.length+1)
										if prevFound == openSize-1 {
											prevFound = i //they will have been swapped
										}
										break searchMatch

									} else {
										openSize, _, _ = align.removeOpenState(prevFound, openSize, resultsSize, s2.length+1)
									}
								}
							}
							found = j
							prevFound = i
							//extend the chain
							ns := align.popState()
							ns.prev = s
							ns.aPos = j
							ns.bPos = bIndex
							ns.aGapIndex = j+2
							ns.aGap = aRed[j+1]
							ns.bGap = bSegments[bIndex+1]
							ns.length = s.length+1
							align.open[i] = ns //replace the chain with the new end
							if (ns.length*2)/3 > minMatches {
								minMatches = (ns.length*2)/3
								maxBIndex = len(bSegments)-minMatches*2+1
							}
							break searchMatch
						}
						g += aRed[j+1] + k
						if debug && g > maxGap {
							log.Println("Stopping at gap index",j+2,"once we reached gap of",g)
						}
					}
				}
				//test for "not enough seeds left"
				if s.length + (len(bSegments) - bIndex) < minMatches {
					openSize, resultsSize, minMatches = align.removeOpenState(i,minMatches, openSize, resultsSize)
				} else {
					//extend the bGap
					s.bGap += bSegments[bIndex+1] + k
				}
			}
		}
		bOffset = 0
		//and look for chains starting here
		if bIndex <= maxBIndex {
			for i := 0; i < initialSize; i++ {
				s := align.initials[i]
				aPos := s.aPos
				if aPos != found && aRed[aPos] == bSeed {
					if found != -1 {
						//do an extra check for an existing state
						for j := 0; j < openSize; j++ {
							if align.open[j].bPos == bIndex && align.open[j].aPos == aPos {
								found = aPos
								break
							}
						}
					}
					if found == aPos || openSize >= len(align.open) {
						continue
					}
					//start a new chain, add to the open list
					ns := align.popState()
					ns.aPos = s.aPos
					ns.bPos = bIndex
					ns.aGapIndex = s.aPos+2
					ns.aGap = aRed[s.aPos+1]
					ns.bGap = bSegments[bIndex+1]
					ns.length = 1
					ns.prev = nil
					if openSize >= len(align.open) {
						for m, op := range align.open[:openSize] {
							log.Println(m,":",op.aPos,op.bPos," len:",op.length)
						}
						log.Println(aRed)
					}
					align.open[openSize] = ns
					openSize++
				}
				//TODO: drop any impossible start positions
			}
		}
		if debug {
			for i := 0; i < openSize; i++ {
				s := align.open[i]
				log.Println("State",i,": ",aMapping[s.aPos/2],"-",s.bPos/2," gaps: ",s.aGap,s.bGap," length:",s.length)
			}
		}
	}
	//add any good ones from the open list
	//log.Println(openSize,"states left at the end.",resultsSize,"already in results.")
	for i := 0; i < openSize; i++ {
		s := align.open[i]
		//log.Println(i,s.aPos,"-",s.bPos,s.length)
		if s.length >= minMatches {
			align.results[resultsSize] = s
			resultsSize++
		}
	}
	if resultsSize == 0 {
		return nil
	}
	matches := make([]*seeds.SeedMatch, 0, resultsSize)
	for i := resultsSize-1; i >= 0; i-- {
		r := align.results[i].extractMatch(aMapping)
		r.SeqA = a
		r.SeqB = b
		matches = append(matches, r)
	}
	return matches
}

func (align *seedAligner) PairwiseBestAlignment(a, b *seeds.SeedSequence, aSet, bSet *util.IntSet, minMatches, k int) *seeds.SeedMatch {
	return nil
}
