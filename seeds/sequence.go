package seeds

import (
	"fmt"
	"github.com/jteutenberg/downpore/util"
	"log"
	"sort"
)

type SeedSequence struct {
	segments          []int //gap size, kmer id interleaved
	id                int
	name              *string
	length            int
	offset            int //offset in bases from the beginning of the original sequence
	inset             int //offset in bases from the end of the original sequence
	reverseComplement *SeedSequence
	rc                bool          //whether this has been reverse-complemented from the original read
	Parent            *SeedSequence //a containing seed-sequence (if applicable)
}

//SeedMatch describes a matching between a subsequence of seeds in two SeedSequences
//It includes a sparse alignment of all exact matching seeds
type SeedMatch struct {
	MatchA                 []int
	MatchB                 []int
	MismatchCount          int
	SeqA                   *SeedSequence
	SeqB                   *SeedSequence
	QueryID                int  //optional
	ReverseComplementQuery bool //whether this match was to the reverse-complement of the original query
}

//LoadSequence produces an otherwise "null" SeedSequence wrapping the giving segments. Useful for testing.
func LoadSequence(segments []int) *SeedSequence {
	return &SeedSequence{segments: segments}
}

//SubSequence creates a new SeedSequence, keeping the given start and end seed indices
//This subsequence shares segment data with its parent
func (s *SeedSequence) SubSequence(start, end, length, offset, inset int) *SeedSequence {
	//e.g. SubSequence 0,0 keeps segments{frontOffset,seed 0, backOffset} <- 3 values
	subs := SeedSequence{segments: s.segments[start*2 : end*2+3], length: length, offset: offset, inset: inset, rc: s.rc, id: s.id, Parent: s}
	return &subs
}

//Trim keeps any seed between and including the specified ones, or those within the offset bases of them
//It returns the number of seeds trimmed from the front
func (s *SeedSequence) Trimmed(startOffset, startSeed, endOffset, endSeed, k int) (*SeedSequence, int) {
	//include any earlier seeds within the range of the bases of startOffset
	for startSeed > 0 && startOffset >= s.segments[startSeed*2]+k {
		startOffset -= s.segments[startSeed*2] + k
		startSeed--
	}
	//and step forward over seeds at the end too
	numSeeds := len(s.segments) / 2
	for endSeed < numSeeds-1 && endOffset >= s.segments[endSeed*2+2]+k {
		endOffset -= s.segments[endSeed*2+2] + k
		endSeed++
	}
	//offset and inset are based on the provided seed + base offsets. At this point start/endOffset are positive.
	offset := s.GetSeedOffset(startSeed, k) - startOffset
	inset := s.GetSeedOffsetFromEnd(endSeed, k) - endOffset
	var trimmed *SeedSequence
	if s.rc {
		//insets and offsets are on opposite ends
		trimmed = s.SubSequence(startSeed, endSeed, s.length-offset-inset, s.offset+inset, s.inset+offset)
	} else {
		trimmed = s.SubSequence(startSeed, endSeed, s.length-offset-inset, s.offset+offset, s.inset+inset)
	}
	segments := make([]int, len(trimmed.segments), len(trimmed.segments))
	copy(segments, trimmed.segments)
	segments[0] = startOffset
	segments[len(trimmed.segments)-1] = endOffset
	trimmed.segments = segments
	return trimmed, startSeed
}

//Reduced creates a version of this seed sequence only containing whitelisted seeds
func (s *SeedSequence) Reduced(whitelist *util.IntSet, k int, minSeeds int, makeIndex bool) (reduced *SeedSequence, index []int) {
	count := 0
	n := len(s.segments)
	prev := -1
	for i := 1; i < n; i += 2 {
		next := s.segments[i]
		if next != prev && whitelist.Contains(uint(next)) {
			count++
			prev = next
		}
	}
	if count < minSeeds {
		return nil, nil
	}
	segs := make([]int, count*2+1)
	offset := s.segments[0]
	if makeIndex {
		index = make([]int, count)
	}
	prev = -1
	j := 0
	for i := 1; i < n; i += 2 {
		seed := s.segments[i]
		if prev != seed && whitelist.Contains(uint(seed)) {
			segs[j] = offset
			segs[j+1] = seed
			if makeIndex {
				index[j/2] = i / 2 //map new seed position to original position
			}
			j += 2
			offset = s.segments[i+1]
			prev = seed
		} else {
			offset += s.segments[i+1] + k
		}
	}
	segs[j] = offset
	return &SeedSequence{segments: segs, length: s.length, offset: s.offset, inset: s.inset, rc: s.rc, id: s.id, Parent: s}, index
}

func ReverseComplement(seed uint, k uint) uint {
	rc := uint(0)
	for j := uint(0); j < k; j++ {
		rc = (rc << 2) | ((seed ^ 3) & 3)
		seed = seed >> 2
	}
	return rc
}

func (s *SeedSequence) ReverseComplement(k int, index *SeedIndex) *SeedSequence {
	if s.reverseComplement != nil {
		return s.reverseComplement
	}
	n := len(s.segments)
	seg := make([]int, n, n)
	n--
	for i, seed := range s.segments {
		if i&1 == 0 {
			seg[n-i] = seed
		} else {
			seed = index.seedMap[seed]
			//A->T is 00 -> 11
			//C->G is 01 -> 10
			rc := 0
			for j := 0; j < k; j++ {
				rc = (rc << 2) | ((seed ^ 3) & 3)
				seed = seed >> 2
			}
			seg[n-i] = index.kmerMap[rc]
		}
	}
	//offset/inset remain as indices in the original forward read
	ns := SeedSequence{seg, s.id, nil, s.length, s.offset, s.inset, s, !s.rc, s.Parent}
	return &ns
}

func (s *SeedSequence) IsReverseComplement() bool {
	return s.rc
}

//Shift adds or removes bases before the first seed, adjusting the offset to compensate
func (s *SeedSequence) Shift(bases int) {
	s.segments[0] += bases
	if s.rc {
		s.inset -= bases
	} else {
		s.offset -= bases
	}
}
func (s *SeedSequence) Extend(bases int) {
	s.segments[len(s.segments)-1] += bases
	if s.rc {
		s.inset -= bases
	} else {
		s.offset -= bases
	}
}

//GetOffset gets the number of bases in the original sequence before the first seed begins
func (s *SeedSequence) GetOffset() int {
	return s.offset
}
func (s *SeedSequence) GetInset() int {
	return s.inset
}

//GetLength gets the length of this seed sequence in bases
func (s *SeedSequence) GetLength() int {
	return s.length
}

//offset gives how far b sequence is ahead of its start seed to be in parity with a
//This will match the two seeds at the start position if they are the same
func (a *SeedSequence) MatchFrom(b *SeedSequence, startA int, startB int, offset int, k int) *SeedMatch {
	m := SeedMatch{MatchA: nil, MatchB: nil, SeqA: a, SeqB: b}
	if startB >= len(b.segments)/2 || startA >= len(a.segments)/2 {
		m.MatchA = make([]int, 0, 0)
		m.MatchB = make([]int, 0, 0)
		return &m
	}
	//how much the offsets of two matching seeds can be, in bases
	maxOffsetRatio := 1.5
	minOffsetRatio := 0.66
	gapLimit := len(a.segments) / 10
	if gapLimit < 5 {
		gapLimit = 5
	}
	matchA := make([]int, 0, len(a.segments)/2-startA)
	matchB := make([]int, 0, len(a.segments)/2-startA)

	minBIndex := startB*2 + 1 //the first seed that might match, probably in the past (by offset)
	maxBIndex := minBIndex + gapLimit*2
	offsetB := -offset //original offset is now applied. Ignore from here on in.
	offsetA := 0
	for i := startA*2 + 1; i < len(a.segments); i += 2 {
		//if whitelist.Contains(uint(a.segments[i])) {
		//determine the limits of how far away the next seed match could be
		minOffset := int(minOffsetRatio * float64(offsetA))
		if minOffset < 0 {
			minOffset = int(float64(offsetA) * maxOffsetRatio) // can be further back
		}
		maxOffset := int(maxOffsetRatio * float64(offsetA))
		if maxOffset < k {
			maxOffset = k //close enough to ignore the gap ratio
		}
		//walk up to the permitted offsets (if not there already)
		for offsetB < minOffset && minBIndex < len(b.segments)-1 {
			offsetB += b.segments[minBIndex+1] + k
			minBIndex += 2
		}
		nextBOffset := offsetB
		for j := minBIndex; j < len(b.segments) && j <= maxBIndex; j += 2 {
			if b.segments[j] == a.segments[i] {
				//exact match. TODO: find exact match with best distance? Dynamic programming instead?
				matchA = append(matchA, i/2)
				matchB = append(matchB, j/2)
				offsetA = 0
				offsetB = b.segments[j+1] + k
				minBIndex = j + 2
				maxBIndex = j + gapLimit*2
				break
			}
			//walk forward until the offsets differ too much
			if nextBOffset < minOffset {
				minBIndex += 2 //later seeds will also not match here
				offsetB += b.segments[j+1] + k
			}
			nextBOffset += b.segments[j+1] + k
			if nextBOffset > maxOffset {
				//no plausable matches before we got too far away
				break
			}
		}
		//}
		offsetA += a.segments[i-1] + k
	}
	//TODO: skipped count
	m.MatchA = matchA
	m.MatchB = matchB
	return &m
}

//MatchTo performs a matching back from the given start points, not matching the final pair of seeds (event if they are the same)
func (a *SeedSequence) MatchTo(b *SeedSequence, startA int, startB int, offset int, k int) *SeedMatch {
	m := SeedMatch{MatchA: nil, MatchB: nil, SeqA: a, SeqB: b}
	if startB <= 0 || startA <= 0 {
		m.MatchA = make([]int, 0, 0)
		m.MatchB = make([]int, 0, 0)
		return &m
	}
	if startB*2-1 >= len(b.segments) {
		startB--
	}
	if startA*2-1 >= len(a.segments) {
		startA--
	}

	//how much the offsets of two matching seeds can be, in bases
	maxOffsetRatio := 1.5
	minOffsetRatio := 0.66
	matchA := make([]int, 0, startA)
	matchB := make([]int, 0, startA)

	maxBIndex := startB*2 - 1                //the first seed that might match
	offsetB := offset + b.segments[startB*2] //back to the match, then one seed more
	offsetA := 0
	//Note: offset gives how far back in the sequence we are looking
	for i := startA*2 - 1; i >= 0; i -= 2 {
		offsetA += a.segments[i+1] + k //the offset to the seed at i
		//if whitelist.Contains(uint(a.segments[i])) {
		//determine the limits of how far away the next seed match could be
		minOffset := int(minOffsetRatio * float64(offsetA))
		if minOffset < 0 {
			minOffset = int(float64(offsetA) * maxOffsetRatio) // can be further back
		}
		maxOffset := int(maxOffsetRatio * float64(offsetA))
		if maxOffset < k {
			maxOffset = k //close enough to ignore the gap ratio
		}
		//walk up to the permitted offsets (if not there already)
		for offsetB < minOffset && maxBIndex > 0 {
			offsetB += b.segments[maxBIndex-1] + k
			maxBIndex -= 2
		}
		nextBOffset := offsetB
		for j := maxBIndex; j >= 0; j -= 2 {
			if b.segments[j] == a.segments[i] {
				//exact match
				matchA = append(matchA, i/2)
				matchB = append(matchB, j/2)
				//offsets are to the preceeding seed now
				if j > 0 {
					offsetA = 0 //a.segments[i-1] + k
					offsetB = b.segments[j-1] + k
				}
				maxBIndex = j - 2
				break
			}
			//walk back until the offsets differ too much. One step at a time.
			if nextBOffset < minOffset {
				maxBIndex -= 2 //later seeds will also not match here
				offsetB += b.segments[j-1] + k
			}
			nextBOffset += b.segments[j-1] + k
			if nextBOffset > maxOffset {
				//no plausable matches before we got too far away
				break
			}
		}
		//}
	}
	if len(matchA) == 0 {
		m.MatchA = matchA
		m.MatchB = matchB
		return &m
	}
	//reverse the matches
	for i := 0; i < len(matchA)/2; i++ {
		end := len(matchA) - i - 1
		t := matchA[i]
		matchA[i] = matchA[end]
		matchA[end] = t
		t = matchB[i]
		matchB[i] = matchB[end]
		matchB[end] = t
	}
	//TODO: skipped count
	m.MatchA = matchA
	m.MatchB = matchB
	return &m
}

/*
func (seq *SeedSequence) SingleMatch(query *SeedSequence, seqSet *util.IntSet, minMatch, k int) *SeedMatch {
	//we'll do a single walk across the sequence, so only reduce the query (which should be the shorter of the two)
	q, qIndex := query.Reduced(seqSet, k, minMatch, true)
	if q == nil {
		return nil
	}
	m := seq.topMatch(q, minMatch, k, false)
	if m == nil {
		return nil
	}
	if qIndex != nil {
		for i, pos := range m.MatchA {
			m.MatchA[i] = qIndex[pos]
		}
	}
	m.SeqA = query
	return m
}*/

func (seq *SeedSequence) Match(query *SeedSequence, querySet *util.IntSet, seqSet *util.IntSet, minMatch, k int) []*SeedMatch {
	s := seq
	q := query
	var qIndex []int
	var sIndex []int
	if querySet != nil {
		s, sIndex = seq.Reduced(querySet, k, minMatch, true)
	}
	if seqSet != nil {
		q, qIndex = query.Reduced(seqSet, k, minMatch, true)
	}
	if s == nil || q == nil {
		return nil
	}
	ms := s.dynamicMatch(q, minMatch, k, false)
	if ms != nil {
		for _, m := range ms {
			//convert back from reduced to original sequences' indices
			if qIndex != nil {
				for i, pos := range m.MatchA {
					m.MatchA[i] = qIndex[pos]
				}
			}
			if sIndex != nil {
				for i, pos := range m.MatchB {
					m.MatchB[i] = sIndex[pos]
				}
			}
			m.SeqA = query
			m.SeqB = seq
		}
	}
	return ms
}

//Full: for a given ref*query position, what length is the best chain to there? Then, take the longest of these.
//Sparse: most ref*query are singular so almost linear set
//Given a matching position, *every* chain to it must use that match <- must simplify things
//Also, at a giving matching position you want to only use the longest match... so...
//... should be able to just walk left->right down the query, building up an interleaved set of chains
func (seq *SeedSequence) dynamicMatch(query *SeedSequence, minMatch, k int, debug bool) []*SeedMatch {
	if minMatch == 0 {
		minMatch = 1
	}
	chainsA := make([][]int, len(query.segments)/2) //chain up to each seed. Early chains will often share arrays with later.
	chainsB := make([][]int, len(query.segments)/2)
	var allGoodChains []*SeedMatch
	for qIndex := 1; qIndex < len(query.segments)-minMatch*2+2; qIndex += 2 {
		if query.segments[qIndex-1] < 0 && qIndex > 1 && query.segments[qIndex+1] < 0 && query.segments[qIndex] == query.segments[qIndex-2] && query.segments[qIndex] == query.segments[qIndex+2] {
			continue //internal to closely spaced repeats. Potential massive repeat region has poor information in this representation
		}
		//1. Consider each matching start position with no existing chain
		querySeedIndex := qIndex / 2
		if chainsA[querySeedIndex] != nil {
			continue
		}
		prevSeed := -1 //to help check for repeats: never start chains on internal repeats, corresponding to the check above
		for i := 1; i < len(seq.segments)-minMatch*2+2; i += 2 {
			nextSeed := seq.segments[i]
			if nextSeed == query.segments[qIndex] && nextSeed != prevSeed && (chainsA[querySeedIndex] == nil || chainsB[querySeedIndex][len(chainsB[querySeedIndex])-1] != i/2) { //matching seeds, and either no chain or a chain on a different path
				//  2. Start a new chain with this match
				if debug {
					log.Println("Starting new chain at query", qIndex/2, "sequence", i/2)
				}
				chainsA[querySeedIndex] = make([]int, 1, (len(query.segments)-qIndex)/2)
				chainsB[querySeedIndex] = make([]int, 1, (len(query.segments)-qIndex)/2)
				chainsA[querySeedIndex][0] = querySeedIndex
				chainsB[querySeedIndex][0] = i / 2
				//  3. Extend the chain forward, setting the chain value at each index as matches are made
				chainA, chainB := extendChain(query, seq, chainsA, chainsB, qIndex, i, k, debug)
				if debug {
					log.Println("Got chain:", chainA, chainB)
				}
				//  4. At the end, if the chain is longest so far, and remaining unchained seeds are fewer, return it
				if len(chainA) >= minMatch {
					if allGoodChains == nil {
						allGoodChains = make([]*SeedMatch, 0, 5)
					}
					nextLength := (len(chainA) * 2) / 3
					if nextLength > minMatch {
						minMatch = nextLength
						//remove any chains shorter than this
						for j := len(allGoodChains) - 1; j >= 0; j-- {
							if len(allGoodChains[j].MatchA) < nextLength {
								allGoodChains[j] = allGoodChains[len(allGoodChains)-1]
								allGoodChains[len(allGoodChains)-1] = nil
								allGoodChains = allGoodChains[:len(allGoodChains)-1]
							}
						}
						//TODO: use base coverage rather than chain length
					}
					allGoodChains = append(allGoodChains, &SeedMatch{chainA, chainB, 0, query, seq, -1, false})
					remaining := 0
					for _, c := range chainsA {
						if c == nil {
							remaining++
						}
					}
					if remaining < len(chainA) {
						if debug {
							log.Println("Only", remaining, "seeds left, so this chain is best!")
						}
						return allGoodChains
					}
				}
			}
			prevSeed = nextSeed
		}
	}
	return allGoodChains
}

/*//Similar to dynamicMatch but avoids memory allocation as much as possible, just returning the best match
//This destroys the query sequence.
func (seq *SeedSequence) topMatch(query *SeedSequence, minMatch, k int, debug bool) *SeedMatch {
	if minMatch == 0 {
		minMatch = 1
	}
	var bestChain []int
	for qIndex := 1; qIndex < len(query.segments)-minMatch*2+2; qIndex += 2 {
		qSeed := query.segments[qIndex]
		if qSeed == -1 || (query.segments[qIndex-1] < 0 && qIndex > 1 && query.segments[qIndex+1] < 0 && qSeed == query.segments[qIndex-2] && qSeed == query.segments[qIndex+2]) {
			continue //internal to closely spaced repeats. Potential massive repeat region has poor information in this representation
		}
		//1. Consider each matching start position with no existing chain
		querySeedIndex := qIndex / 2
		prevSeed := -1 //to help check for repeats: never start chains on internal repeats, corresponding to the check above
		for i := 1; i < len(seq.segments)-minMatch*2+2; i += 2 {
			nextSeed := seq.segments[i]
			if nextSeed == qSeed && nextSeed != prevSeed { //matching seeds
				//  2. Start a new chain with this match
				if debug {
					log.Println("Starting new chain at query", qIndex/2, "sequence", i/2)
				}
				chainsA[querySeedIndex] = make([]int, 1, (len(query.segments)-qIndex)/2)
				chainsB[querySeedIndex] = make([]int, 1, (len(query.segments)-qIndex)/2)
				chainsA[querySeedIndex][0] = querySeedIndex
				chainsB[querySeedIndex][0] = i / 2
				//  3. Extend the chain forward, setting the chain value at each index as matches are made
				chainA, chainB := extendChain(query, seq, chainsA, chainsB, qIndex, i, k, debug)
				if debug {
					log.Println("Got chain:", chainA, chainB)
				}
				//  4. At the end, if the chain is longest so far, and remaining unchained seeds are fewer, return it
				if len(chainA) >= minMatch {
					if allGoodChains == nil {
						allGoodChains = make([]*SeedMatch, 0, 5)
					}
					nextLength := (len(chainA) * 2) / 3
					if nextLength > minMatch {
						minMatch = nextLength
						//remove any chains shorter than this
						for j := len(allGoodChains) - 1; j >= 0; j-- {
							if len(allGoodChains[j].MatchA) < nextLength {
								allGoodChains[j] = allGoodChains[len(allGoodChains)-1]
								allGoodChains[len(allGoodChains)-1] = nil
								allGoodChains = allGoodChains[:len(allGoodChains)-1]
							}
						}
						//TODO: use base coverage rather than chain length
					}

}*/

//Extend a chain forward, setting entries of chainsA and chainsB as we go.
//Any shorter chains encountered are overwritten (this should be an exceptional case)
//a and b index are indices in the segments slices, i.e. 2x the seed index
func extendChain(a, b *SeedSequence, chainsA, chainsB [][]int, aIndex, bIndex, k int, debug bool) ([]int, []int) {
	currentChainA := chainsA[aIndex/2]
	currentChainB := chainsB[aIndex/2]
	offsetA := a.segments[aIndex+1]
	offsetB := b.segments[bIndex+1]
	aIndex += 2
	bIndex += 2
	//NOTE: a/bIndex always points to the segment after current offsetA/B
	for aIndex < len(a.segments) && bIndex < len(b.segments) {
		//find the next match from current positions/offsets
		aSeedIndex := aIndex / 2
		//if b is to match this seed, it must be between:
		var minBOffset int
		var maxBOffset int
		if offsetA < 0 {
			minBOffset = -k
			maxBOffset = 0
		} else {
			minBOffset = (offsetA*2)/3 - k
			maxBOffset = (offsetA*3)/2 + k
		}
		//if the next b is too far, start moving a forward
		for maxBOffset < offsetB {
			offsetA += a.segments[aIndex+1] + k
			aIndex += 2
			if aIndex >= len(a.segments) {
				return currentChainA, currentChainB
			}
			aSeedIndex = aIndex / 2
			minBOffset = (offsetA*2)/3 - k
			maxBOffset = (offsetA*3)/2 + k
		}
		//similarly, walk b forward to the minimum offset
		for offsetB < minBOffset {
			offsetB += b.segments[bIndex+1] + k
			bIndex += 2
			if bIndex >= len(b.segments) {
				return currentChainA, currentChainB
			}
		}
		//save b position and offset for later -- the next match may need to start checking from here
		oldBIndex := bIndex
		oldBOffset := offsetB

		if debug {
			log.Println("Starting scan at", aIndex/2, bIndex/2, "with offsets", offsetA, offsetB, "up to max", maxBOffset)
		}
		//now scan for a match, up to the maximum offset
		matched := false
		seedA := a.segments[aIndex]
		for offsetB <= maxBOffset {
			//if a match, reset the offsets and step forward
			if seedA == b.segments[bIndex] {
				//NOTE: by only considering the first match, this becomes approximate rather than optimal matching
				if debug {
					log.Println("HIT! At", aIndex/2, bIndex/2, "of seed", seedA, "when offsets were", offsetA, offsetB)
				}
				//an existing chain up to here
				if chainsA[aSeedIndex] != nil {
					if debug {
						log.Println("Hit an existing chain at", aSeedIndex)
					}
					//if it's the same match, and part of a longer chain
					if bIndex/2 == chainsB[aSeedIndex][len(chainsB[aSeedIndex])-1] && len(chainsA[aSeedIndex]) > len(currentChainA) {
						return currentChainA, currentChainB //they have a better chain already
					}
					if debug {
						log.Println("Ignoring the chain and moving forward.")
					}
					//just overwrite otherwise. A bit tricky to re-use an existing chain's tail, and should be very rare.
					//we've still got the other chain saved. TODO: consider *not* overwriting it in some cases?
				}
				//set the chain for this seed
				currentChainA = append(currentChainA, aSeedIndex)
				chainsA[aSeedIndex] = currentChainA
				currentChainB = append(currentChainB, bIndex/2)
				chainsB[aSeedIndex] = currentChainB
				offsetA = a.segments[aIndex+1]
				offsetB = b.segments[bIndex+1]
				aIndex += 2
				bIndex += 2
				matched = true
				break
			} else {
				offsetB += b.segments[bIndex+1] + k
				bIndex += 2
				if bIndex >= len(b.segments) {
					break //try the next a
				}
			}
		}
		//otherwise, step a and its offset, and reset b to the previous minimum offset position
		if !matched {
			offsetA += a.segments[aIndex+1] + k
			aIndex += 2
			offsetB = oldBOffset
			bIndex = oldBIndex
		}
	}
	return currentChainA, currentChainB
}

type cluster struct {
	target             *SeedSequence //One of the members or a consensus sequence
	support            []int         //support for each seed: number of components that aligned to it
	targetAnchor       int
	targetAnchorOffset int
	components         []*SeedSequence //match to one other member of the cluster, or the consensus sequence
	alignments         []*SeedMatch    //alignments against the consensus
}

func printClusters(cs []*cluster, width int) {
	maxLeft := 0
	maxRight := 0
	k := 10
	anchorBases := make([]int, len(cs), len(cs))
	for i, c := range cs {
		if c.targetAnchor >= c.target.GetNumSeeds() {
			log.Println("ERROR in cluster", i, ": anchor at ", c.targetAnchor, "/", c.target.GetNumSeeds())
		}
		anchorBases[i] = c.target.GetSeedOffset(c.targetAnchor, k)
		if anchorBases[i] > maxLeft {
			maxLeft = anchorBases[i]
		}
		r := c.target.GetSeedOffset(c.target.GetNumSeeds()-1, k) - anchorBases[i]
		if r > maxRight {
			maxRight = r
		}
	}
	f := float64(width) / float64(maxLeft+maxRight)
	for i, c := range cs {
		left := int(f*float64(maxLeft-anchorBases[i]) + 0.5)
		for j := 0; j < left; j++ {
			fmt.Print(" ")
		}
		seedCount := 0
		baseCount := c.target.GetSeedOffset(0, k)
		hasAnchor := false
		for j := 0; j < c.target.GetNumSeeds(); j++ {
			seedCount++
			hasAnchor = hasAnchor || j == c.targetAnchor
			if float64(baseCount)*f > 1.0 {
				if hasAnchor {
					fmt.Print("X")
				} else if seedCount == 0 {
					fmt.Print(".")
				} else if seedCount == 1 {
					fmt.Print("_")
				} else if seedCount < len(c.components)/2 {
					fmt.Print("-")
				} else {
					fmt.Print("+")
				}
				hasAnchor = false
				seedCount = 0
				baseCount -= int(1.0/f + 0.5)
			}
			baseCount += c.target.GetNextSeedOffset(j, k)
		}
		fmt.Println()
	}
}

func makeCluster(first *SeedSequence, anchor, anchorOffset, capacity int) *cluster {
	newCluster := cluster{target: first, targetAnchor: anchor, targetAnchorOffset: anchorOffset, components: make([]*SeedSequence, 0, capacity), alignments: make([]*SeedMatch, 0, capacity), support: nil}
	newCluster.components = append(newCluster.components, first)
	//dummy alignment for the original sequence to itself
	length := len(first.segments) / 2
	al := SeedMatch{SeqA: first, SeqB: first, MatchA: make([]int, length, length), MatchB: make([]int, length, length)}
	for i := 0; i < length; i++ {
		al.MatchA[i] = i
		al.MatchB[i] = i
	}
	newCluster.alignments = append(newCluster.alignments, &al)
	return &newCluster
}

func (c *cluster) intersects(other *cluster) bool {
	for _, s := range c.components {
		for _, t := range other.components {
			if s == t {
				return true
			}
		}
	}
	return false
}

func (c *cluster) isDistinct(others []*cluster) bool {
	for _, other := range others {
		if other != c && c.intersects(other) {
			return false
		}
	}
	return true
}

func (c *cluster) addSequence(m *SeedMatch, k int) []int {
	c.alignments = append(c.alignments, m)
	target, newIndices := m.Merge(k, 1.0/(float64(len(c.components)+1.0)))
	c.target = target
	c.targetAnchor = newIndices[c.targetAnchor]
	//and store the full match
	c.components = append(c.components, m.SeqB)
	if c.support == nil {
		//create and fill with the 1s and 2s
		c.support = make([]int, len(target.segments)/2, len(target.segments)/2)
		for i := 0; i < len(c.support); i++ {
			c.support[i] = 1
		}
		for _, i := range m.MatchA {
			c.support[newIndices[i]] = 2
		}
	} else {
		//expand and increment support
		oldSupport := c.support
		c.support = make([]int, len(target.segments)/2, len(target.segments)/2)
		for i := 0; i < len(c.support); i++ {
			c.support[i] = 1
		}
		for i, s := range oldSupport {
			c.support[newIndices[i]] = s
		}
		//and anything matched to gets another support
		for _, i := range m.MatchA {
			c.support[newIndices[i]]++
		}
	}
	//update alignments
	for _, a := range c.alignments {
		for i, mat := range a.MatchA { //the index to seeds in the consensus
			a.MatchA[i] = newIndices[mat]
		}
		a.SeqA = target
	}
	return newIndices
}

//rationalise cleans up the consensus for a cluster by removing unsupported seeds
func (c *cluster) rationalise(k int, keepEdges bool) {
	length := 0
	newIndices := make([]int, len(c.support), len(c.support))
	//find first non-1 support entry
	for length < len(c.support) && c.support[length] == 1 {
		newIndices[length] = length
		length++
	}
	start := 0
	offset := 0 //bases passed since last good seed
	if !keepEdges {
		start = length
		if c.targetAnchor < length {
			start = c.targetAnchor
			for i := start; i < length; i++ {
				newIndices[i] = i - start
			}
		}
		offset = -c.target.segments[length*2] //ensures we start with zero gap at front
	}
	//do the same from the end
	end := len(c.support) - 1
	for end > 0 && c.support[end] == 1 {
		end--
	}
	//collapse the centre down
	for index := length; index <= end; index++ {
		offset += c.target.segments[index*2] //add the offset up to here
		seed := c.target.segments[index*2+1]
		isAnchor := index == c.targetAnchor
		if c.support[index] == 1 && !isAnchor {
			//ignore this seed
			offset += k
		} else {
			newIndices[index] = length - start
			c.support[length] = c.support[index]
			seg := length * 2
			c.target.segments[seg] = offset
			c.target.segments[seg+1] = seed
			length++
			offset = 0
		}
	}
	//and write the tail
	if keepEdges {
		for index := end + 1; index < len(c.support); index++ {
			c.support[length] = c.support[index]
			seed := c.target.segments[index*2+1]
			newIndices[index] = length - start
			seg := length * 2
			c.target.segments[seg] = c.target.segments[index*2] + offset
			offset = 0
			c.target.segments[seg+1] = seed
			length++
		}
		c.targetAnchor = newIndices[c.targetAnchor]
		c.target.segments[length*2] = 0 //c.target.segments[len(c.target.segments)-1]
		c.target.segments = c.target.segments[:length*2+1]
		c.support = c.support[:length]
	} else {
		c.target.segments[length*2] = 0
		c.target.segments = c.target.segments[start*2 : length*2+1]
		c.support = c.support[start:length]
		c.targetAnchor = newIndices[c.targetAnchor]
	}
	//collapse alignments down.
	for _, a := range c.alignments {
		index := 0
		for i := 0; i < len(a.MatchA); i++ {
			m := a.MatchA[i]
			if !keepEdges && m < start {
				continue
			}
			if (keepEdges && m < start) || newIndices[m] != 0 { //not deleted
				a.MatchA[index] = newIndices[m]
				a.MatchB[index] = a.MatchB[i]
				index++
			}
		}
		a.MatchA = a.MatchA[:index]
		a.MatchB = a.MatchB[:index]
	}
}

//ReverseComplement replaces both SeqA and SeqB, reverses the match and corrects the indices
func (m *SeedMatch) ReverseComplement(k int, index *SeedIndex) {
	m.SeqA = m.SeqA.ReverseComplement(k, index)
	m.SeqB = m.SeqB.ReverseComplement(k, index)
	end := len(m.MatchA) - 1
	lengthA := len(m.SeqA.segments)/2 - 1
	lengthB := len(m.SeqB.segments)/2 - 1
	//reverse the seed match order so it fits the RC sequences
	for i := 0; i < len(m.MatchA)/2; i++ {
		m.MatchA[i], m.MatchA[end-i] = m.MatchA[end-i], m.MatchA[i]
		m.MatchB[i], m.MatchB[end-i] = m.MatchB[end-i], m.MatchB[i]
	}
	//then change their values so the seed positions are correct. These will be ascending order again
	for i := 0; i < len(m.MatchA); i++ {
		m.MatchA[i] = lengthA - m.MatchA[i]
		m.MatchB[i] = lengthB - m.MatchB[i]
	}
}

func (m *SeedMatch) Validate() bool {
	for i := 0; i < len(m.MatchA); i++ {
		if m.SeqA.segments[m.MatchA[i]*2+1] != m.SeqB.segments[m.MatchB[i]*2+1] {
			log.Println("Invalid alignment at index ", i, ": ", m.SeqA.segments[m.MatchA[i]*2+1], m.SeqB.segments[m.MatchB[i]*2+1])
			log.Println(m.MatchA, "\n", m.MatchB)
			log.Println(m.SeqA, "\n", m.SeqB)
			return false
		}
	}
	return true
}

func (m *SeedMatch) GetBasesCovered(k int) (int, int) {
	countA := len(m.MatchA) * k
	countB := countA
	prevA := m.MatchA[0]
	prevB := m.MatchB[0]
	for i, s := range m.MatchA {
		if i == 0 {
			continue
		}
		d1 := m.SeqA.segments[prevA*2+2]
		d2 := m.SeqB.segments[prevB*2+2]
		for j := prevA + 2; j <= s; j++ {
			d1 += m.SeqA.segments[j*2] + k
		}
		s2 := m.MatchB[i]
		for j := prevB + 2; j <= s2; j++ {
			d2 += m.SeqB.segments[j*2] + k
		}
		if d1 < 0 { //overlap
			countA += d1 //so we subtract the overlap
		}
		if d2 < 0 {
			countB += d2
		}
		prevB = s2
		prevA = s
	}
	return countA, countB
}

//isFullMatch tests whether a mismatch exists at the edge of a sequence in the alignment
func isFullMatch(m *SeedMatch) bool {
	seedLimit := 5 //maximum unmatched seeds at an edge
	gapLimit := 10 //maximum unmatched run of seeds away from an edge
	if len(m.MatchA) < seedLimit {
		return false
	}
	mLen := len(m.MatchA)
	//if the overlap is large (in #seeds), give them a bit more leeway
	aLimit := (m.MatchA[mLen-1] - m.MatchA[0]) / 10
	bLimit := (m.MatchB[mLen-1] - m.MatchB[0]) / 10
	if aLimit < seedLimit {
		aLimit = seedLimit
	}
	if bLimit < seedLimit {
		bLimit = seedLimit
	}
	if len(m.MatchA)/5 > gapLimit {
		gapLimit = len(m.MatchA) / 5
	}
	//if both LHS are substantially unmatched
	if m.MatchA[0] >= aLimit && m.MatchB[0] >= bLimit {
		return false
	}
	//or if both RHS do not match...
	if len(m.SeqA.segments)/2-m.MatchA[mLen-1] >= aLimit && len(m.SeqB.segments)/2-m.MatchB[mLen-1] >= bLimit {
		return false
	}
	//finally, test for big gaps
	prevA := m.MatchA[0]
	prevB := m.MatchB[0]
	for i, a := range m.MatchA {
		b := m.MatchB[i]
		if a-prevA > gapLimit && b-prevB > gapLimit {
			return false
		}
		prevA = a
		prevB = b
	}
	return true
}

type seqSorter struct {
	seqs   []*SeedSequence
	others [][]int
	value  []int
}

func (s *seqSorter) setStarts(anchors []int, anchorOffsets []int, k int) {
	for i, seq := range s.seqs {
		s.value[i] = seq.GetSeedOffset(anchors[i], k) + anchorOffsets[i]
	}
}
func (s *seqSorter) setEnds(anchors []int, anchorOffsets []int, k int) {
	for i, seq := range s.seqs {
		s.value[i] = seq.GetSeedOffsetFromEnd(anchors[i], k) - anchorOffsets[i]
	}
}
func (s *seqSorter) Len() int {
	return len(s.seqs)
}
func (s *seqSorter) Less(i, j int) bool {
	return s.value[i] < s.value[j]
}
func (s *seqSorter) Swap(i, j int) {
	s.seqs[i], s.seqs[j] = s.seqs[j], s.seqs[i]
	s.value[i], s.value[j] = s.value[j], s.value[i]
	for _, other := range s.others {
		other[i], other[j] = other[j], other[i]
	}
}

func Consensus(seqs []*SeedSequence, badness []int, anchors, anchorOffsets []int, k int) []*SeedMatch {
	//1. sort by quality (seeds matching original query)
	sorter := &seqSorter{seqs: seqs, others: [][]int{anchors, anchorOffsets}, value: badness}
	sort.Sort(sorter)

	minMatchLength := 5

	//2. combine
	retry := make([]int, 0, len(seqs)) //early failures can be retried once later
	c := makeCluster(seqs[0], anchors[0], anchorOffsets[0], len(seqs))
	//NOTE: every 5 sequences, seeds with only 1 supporting sequence are removed. I.e. approx 20% support expected.
	for i := 1; i < len(seqs); i++ {
		mf := c.target.MatchFrom(seqs[i], c.targetAnchor, anchors[i], anchorOffsets[i]-c.targetAnchorOffset, k)
		//back match
		var mb *SeedMatch
		if len(mf.MatchA) == 0 {
			mb = c.target.MatchTo(seqs[i], c.targetAnchor, anchors[i], anchorOffsets[i]-c.targetAnchorOffset, k)
		} else { //use a shared seed as the anchor
			mb = c.target.MatchTo(seqs[i], mf.MatchA[0], mf.MatchB[0], 0, k)
		}
		if len(mb.MatchA)+len(mf.MatchA) > minMatchLength { //how to pick a threshold for "good enough" match?
			m := SeedMatch{SeqA: mb.SeqA, SeqB: seqs[i]}
			//create the full length match
			m.MatchA = append(mb.MatchA, mf.MatchA...)
			m.MatchB = append(mb.MatchB, mf.MatchB...)
			//then generate the merged consensus seed sequence
			c.addSequence(&m, k)
			if len(c.components)%5 == 0 { //TODO: more things to test here: number of 2+ seeds
				c.rationalise(k, false)
			}
		} else {
			retry = append(retry, i)
		}
	}
	//3. Retry any that failed initially
	for _, i := range retry {
		mf := c.target.MatchFrom(seqs[i], c.targetAnchor, anchors[i], anchorOffsets[i]-c.targetAnchorOffset, k)
		var mb *SeedMatch
		if len(mf.MatchA) == 0 {
			continue
		}
		mb = c.target.MatchTo(seqs[i], mf.MatchA[0], mf.MatchB[0], 0, k)
		if len(mf.MatchA)+len(mb.MatchA) > minMatchLength {
			m := SeedMatch{SeqA: mb.SeqA, SeqB: seqs[i]}
			m.MatchA = append(mb.MatchA, mf.MatchA...)
			m.MatchB = append(mb.MatchB, mf.MatchB...)
			c.addSequence(&m, k)
		}
		if len(c.components)%5 == 0 { //TODO: more things to test here: number of 2+ seeds
			c.rationalise(k, false)
		}
	}

	//4. Re-align all to the consensus and return the results
	result := make([]*SeedMatch, 0, len(c.components))
	if len(c.components) == 1 {
		return result
	}
	if len(c.components)%5 != 0 {
		c.rationalise(k, true)
	}

	totalSupport := 0 //total appearances of the consensus seeds = reads*hits
	for _, s := range c.support {
		totalSupport += s
	}
	requiredSupport := (totalSupport * 5) / len(c.support) //5 seeds with average (mean) support
	for j, s := range c.components {
		//pick a seed as anchor
		anchorA := c.alignments[j].MatchA[len(c.alignments[j].MatchA)/2]
		anchorB := c.alignments[j].MatchB[len(c.alignments[j].MatchB)/2]
		mf := c.target.MatchFrom(s, anchorA, anchorB, 0, k)
		if len(mf.MatchA) > 0 {
			mb := c.target.MatchTo(s, mf.MatchA[0], mf.MatchB[0], 0, k)
			if len(mb.MatchA)+len(mf.MatchA) > minMatchLength {
				m := SeedMatch{SeqA: c.target, SeqB: s, MatchA: append(mb.MatchA, mf.MatchA...), MatchB: append(mb.MatchB, mf.MatchB...)}
				support := 0
				for _, n := range m.MatchA {
					support += c.support[n]
				}
				if support >= requiredSupport {
					result = append(result, &m)
				}
			}
		}
	}
	/*if debug {
		for _, r := range result {
			support := 0
			for _, n := range r.MatchA {
				support += c.support[n]
			}
			log.Println(r.SeqB.GetName(),support,"\n",r.LongString(k))
		}
	}*/
	//update the length of the consensus
	if len(result) > 0 {
		result[0].SeqA.length = result[0].SeqA.GetSeedOffset(result[0].SeqA.GetNumSeeds(), k)
	}
	return result
}

//Merge will combine two seed sequences with the given alignment, maintaining all seeds
//The returned slice holds a mapping from old seed indices to new ones.
func (m *SeedMatch) Merge(k int, bWeight float64) (*SeedSequence, []int) {
	sa := m.SeqA.segments
	sb := m.SeqB.segments
	newAIndices := make([]int, len(sa)/2, len(sa)/2)
	maxSeg := make([]int, 0, len(sa)+len(sb)-len(m.MatchA))
	//append left from the overlap, assuming a shared timeline
	i := m.MatchA[0]*2 - 1
	j := m.MatchB[0]*2 - 1
	offsetA := sa[i+1] //relative offset
	offsetB := sb[j+1]
	for i > 0 || j > 0 {
		//interleave back to the start
		//TODO: what about weighted distances??
		if (offsetA < offsetB || j <= 0) && i > 0 {
			maxSeg = append(maxSeg, offsetA)
			maxSeg = append(maxSeg, sa[i])
			newAIndices[i/2] = (len(maxSeg) - 1) / 2
			i -= 2
			offsetB -= offsetA + k
			offsetA = sa[i+1]
		} else {
			maxSeg = append(maxSeg, offsetB)
			maxSeg = append(maxSeg, sb[j])
			j -= 2
			offsetA -= offsetB + k
			offsetB = sb[j+1]
		}
	}
	//add the initial offset (zero it)
	maxSeg = append(maxSeg, 0)

	//and reverse it
	for i := 0; i < len(maxSeg)/2; i++ {
		s := maxSeg[i]
		maxSeg[i] = maxSeg[len(maxSeg)-1-i]
		maxSeg[len(maxSeg)-1-i] = s
	}
	n := m.MatchA[0] //reverse index map up to here
	//indices currently say how far *back* the seeds are but are in-order
	size := len(maxSeg) / 2
	for i := 0; i < n; i++ {
		newAIndices[i] = size - 1 - newAIndices[i]
	}

	//merge between the matching seeds, assuming linear movement in both sequences
	for n := 0; n < len(m.MatchA)-1; n++ {
		i := m.MatchA[n]*2 + 1
		j := m.MatchB[n]*2 + 1
		i2 := m.MatchA[n+1]*2 + 1
		j2 := m.MatchB[n+1]*2 + 1
		maxSeg = append(maxSeg, sa[i])
		newAIndices[i/2] = (len(maxSeg) - 1) / 2

		if i+2 == i2 && j+2 == j2 {
			//just add the weighted mean distance
			maxSeg = append(maxSeg, int((1.0-bWeight)*float64(sa[i+1])+bWeight*float64(sb[j+1])+0.5))
			continue
		}
		//interleave offset-seed-offset
		aLength := float64(m.SeqA.getSeedOffsetBetween(i/2, i2/2, k))
		bLength := float64(m.SeqB.getSeedOffsetBetween(j/2, j2/2, k))
		aFactor := 1.0 - bWeight + bWeight*bLength/aLength //adjusted a bit towards b
		bFactor := bWeight + (1.0-bWeight)*aLength/bLength //adjusted more towards a, usually
		if aLength < float64(k) && bLength < float64(k) {
			aFactor = 1.0
			bFactor = 1.0
		}
		offsetA = sa[i+1]
		offsetB = sb[j+1]
		if offsetA >= k {
			offsetA = int(float64(sa[i+1])*aFactor + 0.5)
		}
		if offsetB >= k {
			offsetB = int(float64(sb[j+1])*bFactor + 0.5) //weighted offset in bases
		}
		i += 2
		j += 2
		lastOffset := offsetA
		for i < i2 || j < j2 {
			for (offsetA <= offsetB || j >= j2) && i < i2 {
				maxSeg = append(maxSeg, offsetA)
				maxSeg = append(maxSeg, sa[i])
				offsetB -= offsetA + k
				offsetA = sa[i+1]
				if offsetA >= k {
					offsetA = int(float64(sa[i+1])*aFactor + 0.5)
				}
				newAIndices[i/2] = (len(maxSeg) - 1) / 2
				i += 2
				lastOffset = offsetA
			}
			for (offsetB < offsetA || i >= i2) && j < j2 {
				maxSeg = append(maxSeg, offsetB)
				maxSeg = append(maxSeg, sb[j])
				offsetA -= offsetB + k
				offsetB = sb[j+1]
				if offsetB >= k {
					offsetB = int(float64(sb[j+1])*bFactor + 0.5)
				}
				j += 2
				lastOffset = offsetB
			}
		}
		maxSeg = append(maxSeg, lastOffset)
	}
	//append any unaligned tail
	i = m.MatchA[len(m.MatchA)-1]*2 + 1
	j = m.MatchB[len(m.MatchA)-1]*2 + 1
	maxSeg = append(maxSeg, sa[i])
	newAIndices[i/2] = (len(maxSeg) - 1) / 2
	i += 2
	j += 2
	offsetA = sa[i-1] //relative offset
	offsetB = sb[j-1]
	for i < len(sa) || j < len(sb) {
		if (offsetA < offsetB || j >= len(sb)) && i < len(sa) {
			maxSeg = append(maxSeg, offsetA)
			maxSeg = append(maxSeg, sa[i])
			newAIndices[i/2] = (len(maxSeg) - 1) / 2
			i += 2
			offsetB -= offsetA + k
			offsetA = sa[i-1]
		} else {
			maxSeg = append(maxSeg, offsetB)
			maxSeg = append(maxSeg, sb[j])
			j += 2
			offsetA -= offsetB + k
			offsetB = sb[j-1]
		}
	}
	//add the final offset (zero it)
	maxSeg = append(maxSeg, 0)

	//TODO: length: calculate from the segments
	//id should be a new consensus id? To be set later
	s := SeedSequence{segments: maxSeg, length: 0, id: -1, offset: 0, rc: false}
	return &s, newAIndices
}

//GetBaseIndex gets the position in seed + offset bases in sequence B of a seed found in
//sequence A. The returned index is that before the match, with offset giving the following bases.
//Note that this is not suitable for use with Trim which expects bases *prior* to the match.
//The base offset includes k bases for the returned seed's length
//distance is the distance in bases of sequence b from the last matching seed
func (m *SeedMatch) GetBaseIndex(aIndex int, k int) (index int, bases int, distance int) {
	//find the before and after seed matches..
	before := 0
	for before < len(m.MatchA) && m.MatchA[before] <= aIndex {
		before++
	}
	if before == 0 {
		//special case: looking for a spot before the first matching seed
		offset := 0
		for i := m.MatchA[0]; i > aIndex; i-- {
			offset += m.SeqA.segments[i*2] + k
		}
		//offset bases before index A. Count back from the first match
		bIndex := m.MatchB[0]
		distance = 0
		for i := bIndex * 2; i > 0 && offset > 0; i -= 2 {
			offset -= m.SeqB.segments[i] + k
			distance += m.SeqB.segments[i] + k
			bIndex--
		}
		if bIndex == 0 {
			return 0, -offset, distance + offset //the index in a is before us, so this will probably be at a negative offset
		}
		return bIndex, -offset, distance //-offset are the extra b bases we removed past the beginning of a
	}
	before--
	bIndex := m.MatchB[before]
	//test for exact match
	if aIndex == m.MatchA[before] {
		return bIndex, 0, 0
	}
	//then use the remaining bases from sequence a
	offset := 0
	for i := m.MatchA[before] + 1; i <= aIndex; i++ {
		offset += m.SeqA.segments[i*2] + k //offset up to this later seed
	}
	distance = 0
	//and remove offset for any additional seeds that appear in b
	for i := bIndex*2 + 2; i < len(m.SeqB.segments) && offset >= m.SeqB.segments[i]; i += 2 {
		offset -= m.SeqB.segments[i] + k
		distance += m.SeqB.segments[i] + k
		bIndex++
	}
	if bIndex >= len(m.SeqB.segments)/2 { //we ran over the end
		return bIndex - 1, offset, distance + offset
	}
	return bIndex, offset, distance + offset
}

func (s *SeedSequence) GetSeedOffset(index int, k int) int {
	index = index*2 + 1
	offset := s.segments[0]
	for i := 2; i < index; i += 2 {
		offset += s.segments[i] + k
	}
	return offset
}

//GetSeedAtOffsetFrom is an inverse of the "GetSeedOffset" functions. It finds the farthest seed within
//the given offset from the source seed
func (s *SeedSequence) GetSeedAtOffsetFrom(offset int, index int, k int) int {
	if offset > 0 {
		i := index*2 + 1
		offset -= s.segments[i+1] + k
		for i < len(s.segments)-1 && offset > 0 {
			i += 2
			offset -= s.segments[i+1] + k
		}
		return i
	}
	i := index*2 + 1
	offset += s.segments[i-1] + k
	for i > 1 && offset < 0 {
		i += 2
		offset += s.segments[i+1] + k
	}
	return i
}

func (s *SeedSequence) GetSeedOffsetFromEnd(index int, k int) int {
	index = index*2 + 1
	offset := s.segments[len(s.segments)-1]
	for i := len(s.segments) - 3; i > index; i -= 2 {
		offset += s.segments[i] + k
	}
	return offset
}

func (s *SeedSequence) GetNextSeedOffset(index, k int) int {
	return s.segments[index*2+2] + k
}

func (s *SeedSequence) GetSeed(index int) int {
	return s.segments[index*2+1]
}

func (s *SeedSequence) GetSegments() []int {
	return s.segments //a shame to expose this. Any alternatives?
}

func (s *SeedSequence) getSeedOffsetBetween(indexA, indexB, k int) int {
	index := indexA*2 + 3 //first seed after indexA
	indexB = indexB*2 + 1
	offset := s.segments[index-1]
	for ; index < indexB; index += 2 {
		offset += s.segments[index+1] + k
	}
	return offset
}

//GetAIndices gets the start and end bases in the original sequence
func (m *SeedMatch) GetAIndices(k int) (start int, end int) {
	start = m.SeqA.segments[0] + m.SeqA.offset
	startA := m.MatchA[0]
	endA := m.MatchA[len(m.MatchA)-1]
	for i := 1; i < startA*2+1; i += 2 {
		start += m.SeqA.segments[i+1] + k
	}
	end = start
	for i := startA*2 + 1; i < endA*2+1; i += 2 {
		end += m.SeqA.segments[i-1] + k
	}
	return start, end
}

func (m *SeedMatch) GetBIndices(k int) (start int, end int) {
	start = m.SeqB.segments[0] + m.SeqB.offset
	startB := m.MatchB[0]
	endB := m.MatchB[len(m.MatchB)-1]
	for i := 1; i < startB*2+1; i += 2 {
		start += m.SeqB.segments[i+1] + k
	}
	end = start
	for i := startB*2 + 1; i < endB*2+1; i += 2 {
		end += m.SeqB.segments[i-1] + k
	}
	return start, end
}

func (m *SeedMatch) String() string {
	s1 := ""
	s2 := "\n"
	for i := 0; i < len(m.MatchA); i++ {
		s1 = fmt.Sprint(s1, " ", m.MatchA[i])
		s2 = fmt.Sprint(s2, " ", m.MatchB[i])
	}
	return s1 + s2
}

func (m *SeedMatch) LongString(seedIndex *SeedIndex) string {
	s := ""
	aOff := 0
	bOff := 0
	k := int(seedIndex.GetSeedLength())
	for i := 0; i < len(m.MatchA); i++ {
		pA := aOff + k
		pB := bOff + k
		aOff = m.SeqA.GetSeedOffset(m.MatchA[i], k)
		bOff = m.SeqB.GetSeedOffset(m.MatchB[i], k)
		s = fmt.Sprint(s, " <", aOff-pA, "/", bOff-pB, "> ", seedIndex.SeedString(m.SeqA.segments[m.MatchA[i]*2+1]))
	}
	return s
}

func (seq *SeedSequence) GetID() int {
	return seq.id
}

func (seq *SeedSequence) GetName() string {
	parent := seq
	for parent.Parent != nil {
		parent = parent.Parent
	}
	if parent.name == nil {
		return fmt.Sprint(parent.id)
	}
	return *(parent.name)
}

func (seq *SeedSequence) SetID(id int) {
	seq.id = id
}

func (seq *SeedSequence) Len() int {
	return seq.length
}

func (s *SeedSequence) GetNumSeeds() int {
	return len(s.segments) / 2
}

func (seq *SeedSequence) String() string {
	s := fmt.Sprint(seq.id, ":")
	for i, v := range seq.segments {
		if i&1 == 0 {
			s = s + fmt.Sprint("<", v, ">")
		} else {
			s = s + fmt.Sprint(" ", v, " ")
		}
	}
	return s
}
func (seq *SeedSequence) LongString(seedIndex *SeedIndex) string {
	s := fmt.Sprint(seq.id, ":")
	for i, v := range seq.segments {
		if i&1 == 0 {
			s = s + fmt.Sprint("<", v, ">")
		} else {
			s = s + fmt.Sprint(" ", seedIndex.SeedString(v), " ")
		}
	}
	return s
}
func (seq *SeedSequence) SubString(start int) string {
	s := fmt.Sprint(seq.id, ": ...")
	for i, v := range seq.segments {
		if i/2 < start {
			continue
		}
		if i&1 == 0 {
			s = s + fmt.Sprint("<", v, ">")
		} else {
			s = s + fmt.Sprint(" ", v, " ")
		}
	}
	return s
}
