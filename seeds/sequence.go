package seeds

import (
	"fmt"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/util"
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
		//fmt.Println("shifting end up from ",endSeed," because offset ",endOffset," > ",s.segments[endSeed*2+2])
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

func ReverseComplement(seed uint, k uint) uint {
	rc := uint(0)
	for j := uint(0); j < k; j++ {
		rc = (rc << 2) | ((seed ^ 3) & 3)
		seed = seed >> 2
	}
	return rc
}

func (s *SeedSequence) ReverseComplement(k int) *SeedSequence {
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
			//A->T is 00 -> 11
			//C->G is 01 -> 10
			rc := 0
			for j := 0; j < k; j++ {
				rc = (rc << 2) | ((seed ^ 3) & 3)
				seed = seed >> 2
			}
			seg[n-i] = rc
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
		nextBOffset := offsetB
		//determine the limits of how far away the next seed match could be
		minOffset := int(minOffsetRatio * float64(offsetA))
		if minOffset < 0 {
			minOffset = int(float64(offsetA) * maxOffsetRatio) // can be further back
		}
		maxOffset := int(maxOffsetRatio * float64(offsetA))
		if maxOffset < k {
			maxOffset = k //close enough to ignore the gap ratio
		}
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
		nextBOffset := offsetB
		//determine the limits of how far away the next seed match could be
		minOffset := int(minOffsetRatio * float64(offsetA))
		if minOffset < 0 {
			minOffset = int(float64(offsetA) * maxOffsetRatio) // can be further back
		}
		maxOffset := int(maxOffsetRatio * float64(offsetA))
		if maxOffset < k {
			maxOffset = k //close enough to ignore the gap ratio
		}
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

//Match tests for a match of two gapped-seed sequences, with full-length match on the "query"
//Seeds not in the provided seedSet are ignored
func (seq *SeedSequence) Match(query *SeedSequence, seedSet *util.IntSet, minMatch, k int) *SeedMatch {
	maxOffsetRatio := 1.5
	minOffsetRatio := 0.66
	matchA := make([]int, 0, minMatch)
	matchB := make([]int, 0, minMatch)

	var bestMatch *SeedMatch
	for i := 1; i < len(seq.segments); i += 2 {
		next := seq.segments[i]
		if seedSet != nil && !seedSet.Contains(uint(next)) {
			continue
		}
		qStart := -1
		for j := 1; j < len(query.segments)-minMatch; j += 2 {
			if next == query.segments[j] {
				qStart = j
				break
			}
		}
		if qStart != -1 {
			matchA = append(matchA, qStart/2)
			matchB = append(matchB, i/2)
			//test completion here
			qOffset := 0 //how many bases since last kmer
			offset := seq.segments[i+1]
			minIndex := i + 2
			//step forward through the query looking for matches
			for qIndex := qStart + 2; qIndex < len(query.segments); qIndex += 2 {
				qOffset += query.segments[qIndex-1]
				nextOffset := offset
				minOffset := int(minOffsetRatio * float64(qOffset))
				if minOffset < 0 {
					minOffset = int(float64(qOffset) * maxOffsetRatio) // can be further back
				}
				maxOffset := int(maxOffsetRatio * float64(qOffset))
				if maxOffset < k {
					maxOffset = k //close enough to ignore the gap ratio
				}
				for j := minIndex; j < len(seq.segments) && nextOffset < maxOffset; j += 2 {
					if query.segments[qIndex] == seq.segments[j] {
						matchA = append(matchA, qIndex/2)
						matchB = append(matchB, j/2)
						qOffset = query.segments[qIndex+1] + k
						offset = seq.segments[j+1] + k
						minIndex = j + 2
						break
					}
					//walk forward until the offsets differ too much
					if nextOffset < minOffset {
						minIndex += 2 //later seeds will also not match here
						offset += seq.segments[j+1] + k
					}
					nextOffset += seq.segments[j+1] + k
					if nextOffset > maxOffset {
						//no plausable matches before we got too far away
						break
					}
				}
				//test for an end condition (no chance of matching enough)
				if len(matchA)+(len(query.segments)-qIndex-2)/2 < minMatch {
					break
				}
			}
			if len(matchA) >= minMatch && (bestMatch == nil || len(matchA) > len(bestMatch.MatchA)) {
				skipped := len(query.segments)/2 - len(matchA) //TODO:ignores seq seeds.
				//distance measure where any missed seeds are assumed to be 2-3 base errors, i.e. divide out k
				sm := SeedMatch{matchA, matchB, skipped, query, seq, -1, false}
				if bestMatch != nil {
					matchA = bestMatch.MatchA[:0]
					matchB = bestMatch.MatchB[:0]
				} else {
					matchA = make([]int, 0, minMatch)
					matchB = make([]int, 0, minMatch)
				}
				bestMatch = &sm
			}
			matchA = matchA[:0]
			matchB = matchB[:0]
		}
	}
	return bestMatch
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
			fmt.Println("ERROR in cluster", i, ": anchor at ", c.targetAnchor, "/", c.target.GetNumSeeds())
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
		//fmt.Println("After:",a.MatchA)
		//a.Validate()
	}
}

//ReverseComplement replaces both SeqA and SeqB, reverses the match and corrects the indices
func (m *SeedMatch) ReverseComplement(k int) {
	m.SeqA = m.SeqA.ReverseComplement(k)
	m.SeqB = m.SeqB.ReverseComplement(k)
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
			fmt.Println("Invalid alignment at index ", i, ": ", m.SeqA.segments[m.MatchA[i]*2+1], m.SeqB.segments[m.MatchB[i]*2+1])
			fmt.Println(m.MatchA, "\n", m.MatchB)
			fmt.Println(m.SeqA, "\n", m.SeqB)
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
	//1. sort by quality
	sorter := &seqSorter{seqs: seqs, others: [][]int{anchors, anchorOffsets}, value: badness}
	sort.Sort(sorter)

	//2. combine
	retry := make([]int, 0, len(seqs)) //early failures can be retried once later
	c := makeCluster(seqs[0], anchors[0], anchorOffsets[0], len(seqs))
	for i := 1; i < len(seqs); i++ {
		mf := c.target.MatchFrom(seqs[i], c.targetAnchor, anchors[i], anchorOffsets[i]-c.targetAnchorOffset, k)
		//back match
		var mb *SeedMatch
		if len(mf.MatchA) == 0 {
			mb = c.target.MatchTo(seqs[i], c.targetAnchor, anchors[i], anchorOffsets[i]-c.targetAnchorOffset, k)
		} else { //use a shared seed as the anchor
			mb = c.target.MatchTo(seqs[i], mf.MatchA[0], mf.MatchB[0], 0, k)
		}
		if len(mb.MatchA)+len(mf.MatchA) > 5 { //how to pick a threshold for "good enough" match?
			m := SeedMatch{SeqA: mb.SeqA, SeqB: seqs[i]}
			//create the full length match
			m.MatchA = append(mb.MatchA, mf.MatchA...)
			m.MatchB = append(mb.MatchB, mf.MatchB...)
			//then generate the merged consensus seed sequence
			c.addSequence(&m, k)
			//fmt.Println("After adding",i,"(length",m.SeqB.Len(),") length is",c.target.GetSeedOffset(c.target.GetNumSeeds(),k))
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
		if len(mf.MatchA)+len(mb.MatchA) > 5 {
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
	for j, s := range c.components {
		//pick a seed as anchor
		anchorA := c.alignments[j].MatchA[len(c.alignments[j].MatchA)/2]
		anchorB := c.alignments[j].MatchB[len(c.alignments[j].MatchB)/2]
		mf := c.target.MatchFrom(s, anchorA, anchorB, 0, k)
		mb := c.target.MatchTo(s, mf.MatchA[0], mf.MatchB[0], 0, k)
		//matchLen := len(mf.MatchA) + len(mb.MatchA)
		//fmt.Println("Original match ",len(c.alignments[j].MatchA),"final match",matchLen)
		m := SeedMatch{SeqA: c.target, SeqB: s, MatchA: append(mb.MatchA, mf.MatchA...), MatchB: append(mb.MatchB, mf.MatchB...)}
		result = append(result, &m)
		//fmt.Println(i,": seq",s.id,":",m.LongString(k))
	}
	//update the length of the consensus
	if len(result) > 0 {
		result[0].SeqA.length = result[0].SeqA.GetSeedOffset(result[0].SeqA.GetNumSeeds(), k)
	}
	return result
}

//Cluster performs all-vs-all alignment between the arguments.
//All start/end offsets are relative
func Cluster(leftBases, rightBases int, seqs []*SeedSequence, anchors, anchorOffsets []int, k int) ([][]*SeedMatch, []int) {
	clusters := make([]*cluster, 0, 5)
	// 1. Find approximate start and end orderings (by bases) of each sequence
	sorter := &seqSorter{seqs: seqs, others: [][]int{anchors, anchorOffsets}, value: make([]int, len(seqs), len(seqs))}
	sorter.setStarts(anchors, anchorOffsets, k)
	sort.Sort(sorter)
	startPos := make([]int, len(seqs), len(seqs))
	for i, _ := range seqs {
		startPos[i] = i
	}
	sorter.others = append(sorter.others, startPos)
	sorter.setEnds(anchors, anchorOffsets, k)
	sort.Sort(sorter)
	for i, _ := range seqs {
		sorter.value[i] = -(sorter.others[len(sorter.others)-1][i] + i)
	}
	sort.Sort(sorter) //sorted by "fewest starting before + ending after"
	// 2. Pick the single sequence with minimum start + end offsets (by position in ordering, not base)
	// TODO: 3. Sort remaining sequences by % overlap with chosen one. I.e begin with wholly contained sequences. Break ties by largest first.

	// 4. Perform clustering
	triedRecombine := false
	count := 0
	for i := 0; i < len(seqs); i++ {
		if len(clusters) >= 10 {
			if !triedRecombine {
				triedRecombine = true
				clusters = combineClusters(clusters, seqs, anchors, anchorOffsets, k)
			}
		}
		unmatched := true
		for _, c := range clusters {
			//forward match
			count++
			mf := c.target.MatchFrom(seqs[i], c.targetAnchor, anchors[i], anchorOffsets[i]-c.targetAnchorOffset, k)
			//back match
			var mb *SeedMatch
			if len(mf.MatchA) == 0 {
				mb = c.target.MatchTo(seqs[i], c.targetAnchor, anchors[i], anchorOffsets[i]-c.targetAnchorOffset, k)
			} else { //use a shared seed as the anchor
				mb = c.target.MatchTo(seqs[i], mf.MatchA[0], mf.MatchB[0], 0, k)
			}
			m := SeedMatch{SeqA: mb.SeqA, SeqB: seqs[i]}
			//create the full length match
			m.MatchA = append(mb.MatchA, mf.MatchA...)
			m.MatchB = append(mb.MatchB, mf.MatchB...)
			if isFullMatch(&m) { //if one of this and the consensus matches to the left, and one to the right
				unmatched = false
				//then generate the merged consensus seed sequence
				c.addSequence(&m, k)
				if len(c.components)%5 == 0 {
					c.rationalise(k, true)
				}
			}
		}
		//if no best match, make a new cluster
		if unmatched && len(clusters) < 10 {
			clusters = append(clusters, makeCluster(seqs[i], anchors[i], anchorOffsets[i], len(seqs)-i))
		}
	}
	if len(seqs) > 200 {
		fmt.Println(len(seqs), "sequences made", len(clusters), "clusters using", count, "matchings.")
	}
	clusters = combineClusters(clusters, seqs, anchors, anchorOffsets, k)
	//printClusters(clusters, 100)

	result := make([][]*SeedMatch, 0, len(clusters))
	offsets := make([]int, 0, len(clusters)) //final anchor offsets (in bases)
	for _, c := range clusters {
		//fmt.Println(i,"is distinct:",c.isDistinct(clusters))
		if len(c.components) == 1 || !c.isDistinct(clusters) {
			continue
		}
		//realign each component of a good cluster against its consensus.
		ms := make([]*SeedMatch, 0, len(c.components))
		for j, s := range c.components {
			//pick a seed as anchor
			anchorA := c.alignments[j].MatchA[len(c.alignments[j].MatchA)/2]
			anchorB := c.alignments[j].MatchB[len(c.alignments[j].MatchB)/2]
			//fmt.Println("Match from anchors ",anchorA,"and",anchorB)
			mf := c.target.MatchFrom(s, anchorA, anchorB, 0, k)
			mb := c.target.MatchTo(s, mf.MatchA[0], mf.MatchB[0], 0, k)
			//matchLen := len(mf.MatchA) + len(mb.MatchA)
			//fmt.Println("Original match ",len(c.alignments[j].MatchA),"final match",matchLen)
			m := SeedMatch{SeqA: c.target, SeqB: s, MatchA: append(mb.MatchA, mf.MatchA...), MatchB: append(mb.MatchB, mf.MatchB...)}
			ms = append(ms, &m)
		}
		//find the length of the consensus now
		offset := c.target.GetSeedOffset(c.targetAnchor, k)
		c.target.length = c.target.GetSeedOffsetFromEnd(c.targetAnchor, k) + offset + k
		result = append(result, ms)
		offsets = append(offsets, offset)
	}
	return result, offsets
}

func combineClusters(clusters []*cluster, seqs []*SeedSequence, anchors, anchorOffsets []int, k int) []*cluster {
	for i := 0; i < len(clusters); i++ {
		c := clusters[i]
		if len(c.components) > 1 {
			c.rationalise(k, false)
		}
		//fmt.Println("cluster",i,"has",len(c.target.segments)/2,"seeds from",len(c.components),"sequences.")
		//re-align against earlier clusters.. just in case
		for j := 0; j < i; j++ {
			if len(clusters[j].components) == 1 && len(c.components) == 1 { //this will have been tested on cluster creation
				continue
			}
			//fmt.Println("Comparing against cluster",j,"with anchors",clusters[j].targetAnchor,"and",c.targetAnchor,"out of",len(clusters[j].target.segments)/2,",",len(c.target.segments)/2)
			mf := c.target.MatchFrom(clusters[j].target, c.targetAnchor, clusters[j].targetAnchor, clusters[j].targetAnchorOffset-c.targetAnchorOffset, k)
			mb := c.target.MatchTo(clusters[j].target, c.targetAnchor, clusters[j].targetAnchor, clusters[j].targetAnchorOffset-c.targetAnchorOffset, k)
			m := SeedMatch{SeqA: c.target, SeqB: clusters[j].target}
			//create the full length match
			m.MatchA = append(mb.MatchA, mf.MatchA...)
			m.MatchB = append(mb.MatchB, mf.MatchB...)
			//fmt.Println("Clusters",i,j,"match",len(mf.MatchA),"+", len(mb.MatchA),"against",len(clusters[j].target.segments)/2,"seeds (with",len(c.target.segments)/2,"seeds), from",len(clusters[j].components),"sequences (with my",len(c.components),"sequences). Full = ",isFullMatch(&m), "using anchor",clusters[j].targetAnchor,c.targetAnchor)
			//fmt.Println(m.MatchA)
			//fmt.Println(m.MatchB)
			if isFullMatch(&m) {
				//TODO: merge the shorter one into the larger one
				//printClusters([]*cluster{clusters[j],c},100)
				//each un-shared sequence to the earlier cluster, one by one
				for _, seq := range c.components {
					//check it is not contained in clusters[j] already
					found := false
					for _, otherSeq := range clusters[j].components {
						if otherSeq == seq {
							found = true
							break
						}
					}
					if found {
						continue
					}
					//find a seed shared by the between-consenus alignment and this sequence's alignment to its consensus
					mOffset := 0
					matchedSeedIndex := 0
					targetSeedIndex := clusters[j].targetAnchor
					tOffset := clusters[j].targetAnchorOffset
					for n, s := range seqs {
						if s == seq {
							matchedSeedIndex = anchors[n]
							mOffset = anchorOffsets[n]
							break
						}
					}
					//fmt.Println("Final chosen seeds are:",targetSeedIndex,",",matchedSeedIndex," = ",clusters[j].target.segments[targetSeedIndex*2+1],",",seq.segments[matchedSeedIndex*2+1],"from segments length: ",len(clusters[j].target.segments)/2,",",len(seq.segments)/2)
					mf = clusters[j].target.MatchFrom(seq, targetSeedIndex, matchedSeedIndex, mOffset-tOffset, k)
					if len(mf.MatchA) == 0 {
						//fmt.Println("Match from anchors failed.")
						continue
					}
					//fmt.Println("Match from the seeds: len ",len(mf.MatchA))
					mb = clusters[j].target.MatchTo(seq, mf.MatchA[0], mf.MatchB[0], 0, k)
					newM := SeedMatch{SeqA: clusters[j].target, SeqB: seq}
					//create the full length match
					newM.MatchA = append(mb.MatchA, mf.MatchA...)
					newM.MatchB = append(mb.MatchB, mf.MatchB...)
					// any sequences that fail get discarded (should be rare)
					if isFullMatch(&newM) {
						//fmt.Println("Seq",n,"matched and merged in.")
						//then generate the merged consensus seed sequence
						clusters[j].addSequence(&newM, k)
						if len(c.components)%5 == 0 {
							c.rationalise(k, true)
						}
					}
				}
				// 3. remove this cluster from the list
				clusters[i] = clusters[len(clusters)-1]
				clusters = clusters[:len(clusters)-1]
				i--
				break
			}
		}
	}
	return clusters
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
		//fmt.Println("lengths to next match are",aLength,bLength,"merging with weights",aFactor,bFactor,"=",aFactor+bFactor,"from weighting towards",bWeight)
		//fmt.Println("weighted offsets to next seed are",offsetA,offsetB,"from",sa[i+1],sb[j+1])
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

	/*assert:
	for i := 0; i < len(newAIndices); i++ {
		ind := i*2+1
		newIndex := newAIndices[i]*2+1
		if sa[ind] != maxSeg[newIndex] {
			fmt.Println("BAD REINDEX AT ",i,"/",newAIndices[i],"\n",sa,"\nand\n",sb,"\n",maxSeg,"\nmap is",newAIndices)
		}
	}*/
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
	//fmt.Println("Searching for",aIndex,"seed for its equivalen position in B, k=",k)
	//find the before and after seed matches..
	before := 0
	for before < len(m.MatchA) && m.MatchA[before] <= aIndex {
		before++
	}
	//fmt.Println("First matching seed before it in A is",before," (0= no match before, 1=exact match)")
	if before == 0 {
		//fmt.Println("Prior to first matching seed.")
		//special case: looking for a spot before the first matching seed
		offset := 0
		for i := m.MatchA[0]; i > aIndex; i-- {
			offset += m.SeqA.segments[i*2] + k
		}
		//fmt.Println("Got",offset,"bases between beginnings of",m.MatchA[0],"and",aIndex,"in A")
		//offset bases before index A. Count back from the first match
		bIndex := m.MatchB[0]
		//fmt.Println("Looking about ",offset,"bases back from first match at b index ",bIndex)
		distance = 0
		for i := bIndex * 2; i > 0 && offset > 0; i -= 2 {
			offset -= m.SeqB.segments[i] + k
			distance += m.SeqB.segments[i] + k
			bIndex--
		}
		if bIndex == 0 {
			//fmt.Println("before b's ver first seed! By ",offset," out of our ",m.SeqB.segments[0]," front bases.")
			return 0, -offset, distance + offset //the index in a is before us, so this will probably be at a negative offset
		}
		//fmt.Println("Counting back from first match in B at",m.MatchB[0],"we get to",bIndex,"with",offset,"bases left over")
		return bIndex, -offset, distance //-offset are the extra b bases we removed past the beginning of a
	}
	before--
	bIndex := m.MatchB[before]
	//fmt.Println("The prior matches are at:",m.MatchA[before],",",m.MatchB[before]," in A/B: ",sequence.KmerString(m.SeqA.segments[m.MatchA[before]*2+1],k),"==",sequence.KmerString(m.SeqB.segments[m.MatchB[before]*2+1],k))
	//test for exact match
	if aIndex == m.MatchA[before] {
		return bIndex, 0, 0
	}
	//then use the remaining bases from sequence a
	offset := 0
	for i := m.MatchA[before] + 1; i <= aIndex; i++ {
		offset += m.SeqA.segments[i*2] + k //offset up to this later seed
	}
	//fmt.Println("Counted",offset,"more bases in A to get to",aIndex)
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
	//fmt.Println("B moved up to index",bIndex,"with remaining bases=",offset)
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

func (m *SeedMatch) LongString(k int) string {
	s := ""
	aOff := 0
	bOff := 0
	for i := 0; i < len(m.MatchA); i++ {
		pA := aOff + k
		pB := bOff + k
		aOff = m.SeqA.GetSeedOffset(m.MatchA[i], k)
		bOff = m.SeqB.GetSeedOffset(m.MatchB[i], k)
		s = fmt.Sprint(s, " <", aOff-pA, "/", bOff-pB, "> ", sequence.KmerString(m.SeqA.segments[m.MatchA[i]*2+1], k))
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
func (seq *SeedSequence) LongString(k int) string {
	s := fmt.Sprint(seq.id, ":")
	for i, v := range seq.segments {
		if i&1 == 0 {
			s = s + fmt.Sprint("<", v, ">")
		} else {
			s = s + fmt.Sprint(" ", sequence.KmerString(v, k), " ")
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
