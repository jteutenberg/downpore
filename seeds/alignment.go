package seeds

import (
	//"fmt"
	"github.com/jteutenberg/downpore/util"
	"log"
)

type multiAligner struct {
	segments [][]int
	tempSets []*util.IntSet

	pos []int //seed positions in each sequence
	offset []int //plus bases (0 for exact match)
	gaps []int //distance since last match
}

func NewMultiAligner() *multiAligner {
	return &multiAligner{segments:make([][]int,0,30),tempSets:make([]*util.IntSet,0,30),pos:make([]int,30),offset:make([]int,30),gaps:make([]int,30)}
}

//MSA consensus from approximate shared start position
func (ma *multiAligner) Consensus(seqs []*SeedSequence, k int) (*SeedSequence,[]*SeedMatch) {
	//find all seeds that appear 2+ times
	if len(ma.tempSets) < len(seqs) {
		//find max seed amongst the sequences
		maxSeed := 100
		for _, s := range seqs {
			ms := s.GetMaxSeed()
			if ms > maxSeed {
				maxSeed = ms
			}
		}
		for len(ma.tempSets) < len(seqs) {
			ma.tempSets = append(ma.tempSets, util.NewIntSetCapacity(maxSeed+1))
		}
	}
	for i, set := range ma.tempSets {
		s := seqs[i]
		set.Clear()
		for j := 1; j < len(s.segments); j+=2 {
			set.Add(uint(s.segments[j]))
		}
	}
	useSeeds := util.NewIntSetFromUInts(util.GetSharedIDs(ma.tempSets[:len(seqs)],2,true))

	//reduce the seed sequences to these shared seeds
	seedMap := make([][]int, len(seqs))
	segments := make([][]int, len(seqs))
	for i, s := range seqs {
		red, sm := s.Reduced(useSeeds, k, 1, true)
		if red != nil {
			segments[i] = red.GetSegments()
			seedMap[i] = sm
			//fmt.Println(i,":",segments[i])
		}
	}
	for len(ma.pos) < len(seqs) {
		ma.pos = append(ma.pos,0)
		ma.offset = append(ma.offset,0)
		ma.gaps = append(ma.gaps,0)
	}
	//start state
	for i,_ := range segments {
		ma.pos[i] = -1
		ma.offset[i] = 0
		ma.gaps[i] = 50 //leeway at the start
	}

	consensus := make([]int, 0, useSeeds.Size()*2)
	matches := make([]*SeedMatch, len(seqs))
	for i,_ := range matches {
		if segments[i] != nil {
			matches[i] = &SeedMatch{MatchA:make([]int,0, useSeeds.Size()), MatchB:make([]int,0,useSeeds.Size()),SeqB:seqs[i]}
		}
	}

	supported := make([]int, len(seqs))
	dist := make([]int,len(seqs))
	finished := false
	for !finished {
		fCount := 0
		//get a support count for all next seeds
		near := 100000 //"near" in bases
		for i,segment := range segments {
			p := ma.pos[i]
			supported[i] = 0
			if segment == nil || p >= (len(segment)-1)/2-1 { //e.g. len(segment)==3; if p == 0, too long.
				fCount++
				//already matched all seeds
				continue
			}
			d := segment[p*2+2] - ma.offset[i] //distance to next seed, can actually be negative (slightly behind current position)
			dist[i] = d
			if d < near && d > -k { //don't check ones in the past. They're just here to add support.
				nextSeed := segment[p*2+3]
				minD,maxD := gapRange(d+ma.gaps[i],k)
				minD -= ma.gaps[i]
				maxD -= ma.gaps[i]
				if near > maxD {
					near = maxD
				}
				supported[i] = 1
				for j,segment2 := range segments {
					if segment2 == nil || j == i {
						continue //self-supporting, of course
					}
					p2 := ma.pos[j]+1 //start from next unmatched seed
					//find one more supporter, must be in the future
					if p2 < len(segment2)/2 {
						//actually, this should be checking for an intersection between segment2 range around d (based on its gap) and segment range around d (based on its gap). Perhaps, whichever is higher of the two?
						min2,max2 := gapRange(d+ma.gaps[j],k)
						if min2 > minD {
							min2 = minD
						}
						if max2 < maxD {
							max2 = maxD
						}
						otherD := segment2[p2*2]-ma.offset[j] //remaining distance from current pos to next seed
						for otherD < min2 && p2 < len(segment2)/2 {
							p2++
							otherD += segment2[p2*2]+k //move forward one more seed, adding to the distance
						}
						for otherD < max2 && p2 < len(segment2)/2 { //TODO: missing last seed?
							if segment2[p2*2+1] == nextSeed {
								supported[i]++
								dist[i] += otherD
								break
							}
							p2++
							otherD += segment2[p2*2]+k
						}
					}
				}
			}
		}
		if fCount >= len(seqs) {
			break
		}
		//select minimum distance supported option 
		minseed := -1
		mindist := 0
		minsup := 0
		var minD int
		var maxD int
		for i,d := range dist {
			if supported[i] > 1 {
				d = d/supported[i]
				seed := segments[i][ma.pos[i]*2+3]
				if minseed == -1 || (minseed == seed && supported[i] > minsup) || (minseed != seed && mindist > d) {
					minsup = supported[i]
					mindist = d
					minseed = segments[i][ma.pos[i]*2+3]
					minD,maxD = gapRange(d+ma.gaps[i],k)
					minD -= ma.gaps[i]
					maxD -= ma.gaps[i]
				}
			}
		}
		//fmt.Println(supported,ma.pos[:len(seqs)],ma.offset[:len(seqs)])
		//fmt.Println(dist,"chose",minseed,"distance range",minD,maxD)
		if minseed == -1 {
			//no supports here. Step the shortest gap and continue.
			minIndex := -1
			minDist := 100000
			for i,d := range dist {
				if supported[i] > 1 {
					d = d / supported[i]
				}
				if segments[i] != nil && ma.pos[i] < len(segments)/2 && d < minDist {
					minDist = d
					minIndex = i
				}
			}
			if minIndex == -1 {
				break //out of options
			}
			//increment all by minDist
			for i, segment := range segments {
				if segment != nil {
					ma.gaps[i] += minDist
					ma.offset[i] += minDist
				}
			}
			ma.gaps[minIndex] = 0
			ma.offset[minIndex] = 0
			ma.pos[minIndex]++
			continue
		}
		//add
		consensus = append(consensus, mindist)
		consensus = append(consensus, minseed)

		//build the matchings and step past
		fCount = 0
		for i, segment := range segments {
			if segment == nil {
				fCount++
				continue
			}
			//find the match again
			matchDex := ma.pos[i]+1
			if matchDex < len(segment)/2 {
				//check for expanding min/max bounds
				min2,max2 := gapRange(mindist+ma.gaps[i],k)
				if min2 > minD {
					min2 = minD
				}
				if max2 < maxD {
					max2 = maxD
				}
				otherD := segment[matchDex*2]-ma.offset[i]
				for otherD < min2 && matchDex < len(segment)/2 {
					matchDex++
					otherD += segment[matchDex*2]+k
				}
				found := false
				//if found, add matching and set offset/gap to zero
				for otherD < max2 && matchDex < len(segment)/2 { //TODO: missing last seed?
					if segment[matchDex*2+1] == minseed {
						ma.pos[i] = matchDex
						ma.offset[i] = 0
						ma.gaps[i] = 0
						matches[i].MatchA = append(matches[i].MatchA,len(consensus)/2-1)
						matches[i].MatchB = append(matches[i].MatchB,seedMap[i][matchDex])
						found = true
						break
					}
					matchDex++
					otherD += segment[matchDex*2]+k
				}
				//otherwise, move forward by the distance bases
				if !found {
					ma.gaps[i] += mindist
					ma.offset[i] += mindist
					p := ma.pos[i]
					//step past some seeds if they're well behind
					for p < len(segment)/2 && ma.offset[i] > segment[p*2+2]+50 { //50 bases is a fair buffer
						ma.offset[i] -= segment[p*2+2]+k
						p++
						ma.pos[i]++
					}
					if p >= len(segment)/2 {
						fCount++
					}
				}
			} else {
				fCount++
			}
		}

		//fmt.Println(ma.pos[:len(seqs)],ma.offset[:len(seqs)])

		finished = fCount >= len(seqs)
	}
	consensus = append(consensus,0) //TODO: set to the gap for each?
	seedCons := LoadSequence(consensus)
	for i := len(matches)-1; i >= 0; i-- {
		m := matches[i]
		if m == nil || len(m.MatchA) < 3 {
			matches[i] = matches[len(matches)-1]
			matches = matches[:len(matches)-1]
		} else {
			m.SeqA = seedCons
		}
	}
	return seedCons,matches
}

type Aligner interface {
	PairwiseAlignments(*SeedSequence, *SeedSequence, *util.IntSet, *util.IntSet, int,int,bool) []*SeedMatch
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

func (s *pairState) extractMatch(aMapping []int) *SeedMatch {
	ma := make([]int,s.length)
	mb := make([]int,s.length)
	for s != nil {
		ma[s.length-1] = aMapping[s.aPos/2]
		mb[s.length-1] = s.bPos/2
		s = s.prev
	}
	return &SeedMatch{MatchA:ma, MatchB:mb}
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

func gapRange(gap,k int) (int,int) {
	minGap := (gap*2)/3 - k
	maxGap := (gap*3)/2 + k +1
	if minGap < 0 {
		minGap = -k
		if maxGap < 0 {
			maxGap = 0
		}
	} else if maxGap < 20 {
		maxGap = 20
		minGap = 0
	}
	return minGap,maxGap
}

func (align *seedAligner) PairwiseAlignments(a, b *SeedSequence, aSet, bSet *util.IntSet, minMatches, k int, debug bool) []*SeedMatch {
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
			minGap,maxGap := gapRange(s.bGap,k)
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
	matches := make([]*SeedMatch, 0, resultsSize)
	for i := resultsSize-1; i >= 0; i-- {
		r := align.results[i].extractMatch(aMapping)
		r.SeqA = a
		r.SeqB = b
		matches = append(matches, r)
	}
	return matches
}

