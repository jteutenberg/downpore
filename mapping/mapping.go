package mapping

import (
	"fmt"
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/util"
	"log"
	"sort"
)

type Mapping struct {
	Query       sequence.Sequence
	Start       int
	End         int
	QueryOffset int
	QueryInset  int
	RC          bool
	match       *seeds.SeedMatch
	ids         int //absolute count of identity matches
}

type Mapper interface {
	Map(sequence.Sequence) []*Mapping
	MapWorker(<-chan sequence.Sequence, chan<- []*Mapping, chan<- bool)
	AsString(*Mapping) string
}

type mapper struct {
	index     *seeds.SeedIndex
	reference sequence.Sequence
	edgeSize  int //number of bases to use in queries
	circular  bool
}

type mappingsByPos []*Mapping
type mappingsByQPos []*Mapping
type mappingsByIds []*Mapping

func (m mappingsByPos) Len() int {
	return len(m)
}
func (m mappingsByPos) Less(i, j int) bool {
	return m[i].Start < m[j].Start
}
func (m mappingsByPos) Swap(i, j int) {
	m[i], m[j] = m[j], m[i]
}
func (m mappingsByQPos) Len() int {
	return len(m)
}
func (m mappingsByQPos) Less(i, j int) bool {
	return m[i].QueryOffset < m[j].QueryOffset
}
func (m mappingsByQPos) Swap(i, j int) {
	m[i], m[j] = m[j], m[i]
}
func (m mappingsByIds) Len() int {
	return len(m)
}
func (m mappingsByIds) Less(i, j int) bool {
	return m[i].ids > m[j].ids
}
func (m mappingsByIds) Swap(i, j int) {
	m[i], m[j] = m[j], m[i]
}

func NewMapper(reference sequence.Sequence, circular bool, k uint, kmerValues []float64, seedRate int, edgeSize int, chunkSize int, numWorkers int) Mapper {
	m := mapper{index: seeds.NewSeedIndex(k), reference: reference, edgeSize: edgeSize, circular: circular}
	//walk the sequence adding 1 seed per seedRate bases
	m.index.AddSingleSeeds(reference, seedRate, kmerValues)
	log.Println("Index has", m.index.Size(), "seeds")
	//index chunks of the reference
	chunks := make(chan sequence.Sequence, numWorkers*2)
	seedSeqs := make(chan *seeds.SeedSequence, numWorkers*2)
	done := make(chan bool, numWorkers)
	for i := 0; i < numWorkers; i++ {
		go seeds.AddSequenceWorker(chunks, m.index, seedSeqs, done)
	}
	ind := 0
	go func() {
		for j := 0; j < 10; j++ {
			start := j * chunkSize
			step := chunkSize*10 - edgeSize
			for i := start; i < reference.Len()-chunkSize/2; i += step {
				end := i + chunkSize
				if i >= reference.Len() {
					end = reference.Len()
				}
				//TODO: handle the short sequence at the end?
				chunks <- reference.SubSequence(i, end)
			}

		}
		if circular {
			chunks <- reference.SubSequence(reference.Len()-edgeSize, reference.Len()).Append(0, reference.SubSequence(0, edgeSize), nil)
		}
		close(chunks)
		for i := 0; i < numWorkers; i++ {
			<-done
		}
		close(seedSeqs)
	}()
	for seq := range seedSeqs {
		seq.SetID(ind)
		m.index.AddSequence(seq)
		ind++
	}
	m.index.IndexSequences()
	return &m
}

//AsString prepares a PAF line for this mapping
func (m *mapper) AsString(mapping *Mapping) string {
	rc := "+"
	if mapping.RC {
		rc = "-"
	}
	mappedLength := mapping.End - mapping.Start
	if m.circular && mappedLength < 0 {
		mappedLength = m.reference.Len() - mapping.Start + mapping.End
	}
	return fmt.Sprintf("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t255", mapping.Query.GetName(), mapping.Query.Len(), mapping.QueryOffset, mapping.Query.Len()-mapping.QueryInset, rc, m.reference.GetName(), m.reference.Len(), mapping.Start, mapping.End, mapping.ids, mapping.End-mapping.Start)
}

func updateQuery(ms []*Mapping, q sequence.Sequence) {
	for _, m := range ms {
		m.Query = q
	}
}

//Determines whether two mappings (using the same query) map to the same part of the reference
func (m *mapper) isConsistent(left, right *Mapping) bool {
	if left.RC != right.RC {
		return false
	}
	//distance in query
	expectedDistance := right.QueryOffset - left.Query.Len() + left.QueryInset
	//distance in reference
	var distance int
	if !left.RC {
		distance = right.Start - left.End // - qInset + inset - qOffset + offset //TODO: <- problem: inset/offset change over time
	} else {
		distance = left.Start - right.End // - qInset + inset - qOffset + offset
	}
	if m.circular && distance < -50 {
		distance += m.reference.Len()
	}
	//Distance function: adjacent (say, within 50 bases) is always fine; nearby (say, up to 500 bases) can be 66-150%; distant (say, above 5000 bases) can be 90-111%
	if distance < 50 && expectedDistance < 50 && distance > -50 {
		return true
	}
	if distance < 500 {
		return (expectedDistance < (distance*3)/2 && expectedDistance > (distance*2)/3)
	}
	if distance > 5000 {
		return (expectedDistance < (distance*10)/9 && expectedDistance > (distance*9)/10)
	}
	ratio := float64(distance-500) / 4500.0
	ratio = 3.0/2.0 + ratio*(10.0/9.0-3.0/2.0)
	return distance < int(float64(expectedDistance)*ratio) && distance > int(float64(expectedDistance)/ratio)
}

//maps the (potentially noisy) ends
//if any map end-to-end then these are merged and all others discarded
func (m *mapper) mapEnds(query sequence.Sequence) (openA, openB, matching []*Mapping) {
	openA = m.performMapping(query.SubSequence(0, m.edgeSize))
	openB = m.performMapping(query.SubSequence(query.Len()-m.edgeSize, query.Len()))
	openA = removeDominated(openA, openA, query.Len())
	openB = removeDominated(openB, openB, query.Len())
	updateQuery(openA, query)
	updateQuery(openB, query)
	return m.matchPairs(openA, openB)
}

func (m *mapper) matchPairs(openA, openB []*Mapping) (remainingA, remainingB, matched []*Mapping) {
	for i := len(openA) - 1; i >= 0; i-- {
		ra := openA[i]
		for j := len(openB) - 1; j >= 0; j-- {
			rb := openB[j]
			if m.isConsistent(ra, rb) {
				qOffset := ra.QueryOffset
				qInset := rb.QueryInset
				if ra.RC {
					ra, rb = rb, ra
				}
				combined := &Mapping{Start: ra.Start, End: rb.End, Query: ra.Query, QueryOffset: qOffset, QueryInset: qInset, RC: ra.RC, ids: ra.ids + rb.ids}
				if matched == nil {
					matched = make([]*Mapping, 0, i+1)
				}
				matched = append(matched, combined)
				if ra.RC { //swap back
					ra, rb = rb, ra
				}
				//and remove from the open lists
				openA[i] = openA[len(openA)-1]
				openA = openA[:len(openA)-1]
				openB[j] = openB[len(openB)-1]
				openB = openB[:len(openB)-1]
				break //move to next item of openA
			}
		}
	}
	return openA, openB, matched
}

//search for chimeric boundary given no end-to-end matches
//this assumes any internal sequence matches either some on the left or some on the right, but not both.
func (m *mapper) findSplitPoint(query sequence.Sequence, openA, openB []*Mapping, left, right int) {
	for right-left >= m.edgeSize {
		start := (right + left - m.edgeSize) / 2
		end := start + m.edgeSize
		mid := m.performMapping(query.SubSequence(start, end))
		newLeft := left
		newRight := right
		afterA := 0 //how much of mid got matched with something from openA
		afterB := 0
		for _, mm := range mid {
			mm.Query = query
			for _, ma := range openA {
				if m.isConsistent(ma, mm) {
					ma.QueryInset = mm.QueryInset
					ma.ids += mm.ids
					if ma.RC {
						ma.Start = mm.Start
					} else {
						ma.End = mm.End
					}
					midMatched := query.Len() - mm.QueryInset - mm.QueryOffset
					if midMatched > afterA {
						afterA = midMatched
					}
					if query.Len()-mm.QueryInset > newLeft {
						newLeft = query.Len() - mm.QueryInset
					}
					break
				}
			}
			if afterA < (m.edgeSize*2)/3 {
				//possibility of matching something on the other side too
				for _, mb := range openB {
					if m.isConsistent(mm, mb) {
						mb.QueryOffset = mm.QueryOffset
						mb.ids += mm.ids
						if mb.RC {
							mb.End = mm.End
						} else {
							mb.Start = mm.Start
						}
						midMatched := query.Len() - mm.QueryInset - mm.QueryOffset
						if midMatched > afterB {
							afterB = midMatched
						}
						if mm.QueryOffset < newRight {
							newRight = mm.QueryOffset
						}
						break
					}
				}
			}
		}
		//now see which side was matched (if any, or both!)
		if afterA > 0 && afterB > 0 {
			//this is unusual. Possibly matching in different ways. Probably just found the true mid-point.
			//assume that this is the centre, so we'll just bring the left and right way in, and recurse on them (once).
			empty := make([]*Mapping, 0)
			if newLeft-left > m.edgeSize*2 { //enough left to look at..
				m.findSplitPoint(query, openA, empty, newLeft-m.edgeSize*2, newLeft-m.edgeSize)
			}
			if right-newRight > m.edgeSize*2 {
				m.findSplitPoint(query, empty, openB, newRight+m.edgeSize, newRight+m.edgeSize*2)
			}
			return //and return whatever we found
		}
		if afterA == 0 && afterB == 0 {
			//no match at all?! Possibly junk in the middle. Recurse on each side independently.
			empty := make([]*Mapping, 0)
			if len(openA) > 0 {
				m.findSplitPoint(query, openA, empty, left, start)
			}
			if len(openB) > 0 {
				m.findSplitPoint(query, empty, openB, end, right)
			}
			return
		}
		//otherwise one side matched. Let's continue.
		left = newLeft
		right = newRight
	}
}

func (m *mapper) GetRepeats(start, end int) (starts, ends []int) {
	query := m.reference.SubSequence(start, end)
	hits := m.performMapping(query)
	hits = removeDominated(hits, hits, query.Len())
	starts = make([]int, len(hits))
	ends = make([]int, len(hits))
	for i, mapping := range hits {
		starts[i] = mapping.Start
		ends[i] = mapping.End
	}
	return starts, ends
}

//count A/B are how many of the open list are at an edge (starting at index 0)
//Take further edgeSize steps in, testing for matching pairs and removing redundant hits
func (m *mapper) mapNext(query sequence.Sequence, openA, openB []*Mapping) (newA, newB, matched []*Mapping) {
	var extended []*Mapping

	// Short sequences

	if query.Len() < m.edgeSize*4 {
		//special handling of "short" sequences. Hold on to the edges a bit longer.
		newA = m.performMapping(query.SubSequence(m.edgeSize, query.Len()-m.edgeSize))
		newA = removeDominated(newA, newA, query.Len())
		updateQuery(newA, query)
		openA, newA, extended = m.matchPairs(openA, newA)
		if extended != nil {
			openA = append(newA, extended...) //discard the early-only hits
		} else {
			openA = append(openA, newA...) //keep the early-onlies since nothing looks good..
		}
		//match with B.
		newA, newB, matched = m.matchPairs(openA, openB)
		if matched == nil {
			return newA, newB, matched //return the bits and pieces we found
		}
		return newA[:0], newB[:0], matched //just the good full length one
	}

	// Long sequences

	// 1. Test one in from the edges, see if anything matches up.
	newA = m.performMapping(query.SubSequence(m.edgeSize, m.edgeSize*2))
	newA = removeDominated(newA, newA, query.Len())
	updateQuery(newA, query)
	openA, newA, extended = m.matchPairs(openA, newA)
	openA = append(openA, newA...)
	if extended != nil {
		openA = append(openA, extended...)
	}
	//do the same from the other end
	newB = m.performMapping(query.SubSequence(query.Len()-m.edgeSize*2, query.Len()-m.edgeSize))
	newB = removeDominated(newB, newB, query.Len())
	updateQuery(newB, query)
	openB, newB, extended = m.matchPairs(newB, openB)
	openB = append(openB, newB...)
	if extended != nil {
		openB = append(openB, extended...)
	}

	// see if anything matches full length now
	newA, newB, matched = m.matchPairs(openA, openB)

	// 2. Second (last!) chance: take two steps in from each edge and try again
	if matched == nil {
		//do a second round, just in case we got unlucky with read errors or reference seeds
		if query.Len() > m.edgeSize*5 {
			openA = m.performMapping(query.SubSequence(m.edgeSize*2, m.edgeSize*3))
			openA = removeDominated(openA, openA, query.Len())
			updateQuery(openA, query)
			openA, newA, extended = m.matchPairs(newA, openA)
			if extended != nil {
				openA = append(openA, extended...)
			}
			openA = append(openA, newA...) //keep everything for now. One last chance to match with B
		}
		if query.Len() > m.edgeSize*6 { //room for the other side to extend too
			openB = m.performMapping(query.SubSequence(query.Len()-m.edgeSize*3, query.Len()-m.edgeSize*2))
			openB = removeDominated(openB, openB, query.Len())
			updateQuery(openB, query)
			openB, newB, extended = m.matchPairs(openB, newB)
			if extended != nil {
				openB = append(openB, extended...)
			}
			openB = append(openB, newB...)
		} else {
			openB = newB
		}
		if query.Len() > m.edgeSize*5 {
			newA, newB, matched = m.matchPairs(openA, openB)
		}
	}
	return newA, newB, matched
}

//removeDominated modifies the ordering of open, returning a new slice containing only those mappings not dominated by a member of extended.
//The extended slice must already be ordered by query position
func removeDominated(open, extended []*Mapping, queryLen int) []*Mapping {
	if len(open) == 0 || len(extended) == 0 {
		return open
	}
	sort.Sort(mappingsByQPos(open))
	j := 0 //first relevant index in extended
	//90% contained, 20% less identity match
	toRemove := make([]bool, len(open))
	for i, next := range open {
		for j < len(extended) && queryLen-extended[j].QueryInset < next.QueryOffset {
			j++
		}
		if j == len(extended) {
			return open
		}
		//test all relevant in extended
		dominated := false
		for k := j; !dominated && k < len(extended) && extended[k].QueryOffset < queryLen-next.QueryInset; k++ {
			if extended[k].ids*4 > next.ids*5 {
				//find the bounds of the overlap
				start := next.QueryOffset
				if extended[k].QueryOffset > start {
					start = extended[k].QueryOffset
				}
				end := queryLen - next.QueryInset
				if extended[k].QueryInset > next.QueryInset {
					end = queryLen - extended[k].QueryInset
				}
				dominated = ((end-start)*10 > (queryLen-next.QueryOffset-next.QueryInset)*9)
			}
		}
		toRemove[i] = dominated
	}
	last := len(open) - 1
	for i := last; i >= 0; i-- {
		if toRemove[i] {
			open[i] = open[last]
			last--
		}
	}
	return open[:last+1]
}

func (m *mapper) Map(query sequence.Sequence) []*Mapping {
	var results []*Mapping
	//For short, return the front of the matches (ignoring substantially worse ones)
	if query.Len() <= m.edgeSize*2 {
		results = m.performMapping(query)
		results = removeDominated(results, results, query.Len())
		updateQuery(results, query)
	} else {
		//1. Map the (possibly unreliable) ends of the query
		openA, openB, matched := m.mapEnds(query)

		if matched != nil {
			//the end-to-end mappings are the only mappings of interest
			results = matched
		} else if query.Len() < m.edgeSize*3 { //no room for further checks
			results = append(openA, openB...)
		} else {
			//2. Try two more steps in from each end
			openA, openB, matched = m.mapNext(query, openA, openB)
			if matched != nil {
				//the end-to-end mappings are still the only mappings of interest
				results = matched
			} else {
				//and search for a split point in the remainder, assuming some kind of chimeric situation
				left := m.edgeSize * 2
				right := query.Len() - m.edgeSize*2
				for _, a := range openA {
					if a.QueryInset > left {
						left = a.QueryInset
					}
				}
				left = query.Len() - right
				for _, b := range openB {
					if b.QueryOffset < right {
						right = b.QueryOffset
					}
				}
				m.findSplitPoint(query, openA, openB, left, right)
				//and remove all unpaired ends
				size := query.Len() - m.edgeSize
				for i := len(openA) - 1; i >= 0; i-- {
					if openA[i].QueryInset >= size {
						openA[i] = openA[len(openA)-1]
						openA = openA[:len(openA)-1]
					}
				}
				for i := len(openB) - 1; i >= 0; i-- {
					if openB[i].QueryOffset >= size {
						openB[i] = openB[len(openB)-1]
						openB = openB[:len(openB)-1]
					}
				}
				results = append(openA, openB...)
			}
		}
	}
	return results
}

func (m *mapper) performMapping(query sequence.Sequence) []*Mapping {
	k := int(m.index.GetSeedLength())
	seedQuery := m.index.NewSeedSequence(query)
	rcQuery := m.index.NewSeedSequence(query.ReverseComplement())

	minMatches := seedQuery.GetNumSeeds() / 5
	minRCMatches := rcQuery.GetNumSeeds() / 5
	if minMatches < 5 { //this will be used as a threshold on *unique* seeds (as a set)
		minMatches = 5
	}
	if minRCMatches < 5 {
		minRCMatches = 5
	}
	matchingIndices := m.index.Matches(seedQuery, 0.25) //this threshold counts repeats of seeds in the query
	matchingRCIndices := m.index.Matches(rcQuery, 0.25)
	results := make([]*Mapping, 0, len(matchingIndices)+len(matchingRCIndices))

	//find the best hit
	maxSeed := 0
	for i := 0; i < seedQuery.GetNumSeeds(); i++ {
		s := seedQuery.GetSeed(i)
		if s > maxSeed {
			maxSeed = s
		}
	}
	seedSet := util.NewIntSetCapacity(maxSeed + 1)
	for i := 0; i < seedQuery.GetNumSeeds(); i++ {
		seedSet.Add(uint(seedQuery.GetSeed(i)))
	}
	for _, index := range matchingIndices {
		//1. Does this have enough matching seeds (compared to best so far)?
		matchSet := m.index.GetSeedSet(index)
		if matchSet.CountIntersectionTo(seedSet, minMatches) < uint(minMatches) {
			continue
		}
		//2. Match based on shared seeds
		match := m.index.GetSeedSequence(index)
		seedMatches := match.Match(seedQuery, seedSet, matchSet, minMatches, k)
		if seedMatches != nil {
			for _, seedMatch := range seedMatches {
				start := match.GetOffset() + match.GetSeedOffset(seedMatch.MatchB[0], k)
				end := m.reference.Len() - match.GetInset() - match.GetSeedOffsetFromEnd(seedMatch.MatchB[len(seedMatch.MatchB)-1], k)
				if m.circular && start > m.reference.Len() {
					start -= m.reference.Len() //in the RHS of the circular join
				}
				qOffset := seedQuery.GetSeedOffset(seedMatch.MatchA[0], k)
				qInset := seedQuery.GetSeedOffsetFromEnd(seedMatch.MatchA[len(seedMatch.MatchA)-1], k)
				if qOffset+qInset > (seedQuery.Len()*2)/3 {
					continue
				}
				qOffset += seedQuery.GetOffset()
				qInset += seedQuery.GetInset()
				_, ids := seedMatch.GetBasesCovered(k) //bases matched in the reference
				results = append(results, &Mapping{Start: start, End: end, QueryOffset: qOffset, QueryInset: qInset, RC: false, match: seedMatch, ids: ids})
				limit := (len(seedMatch.MatchA) * 4) / 5
				if limit > minMatches {
					minMatches = limit
				}
				if limit > minRCMatches {
					minRCMatches = limit
				}
			}
		}
	}
	//or reverse-complement match
	seedSet.Clear()
	for i := 0; i < rcQuery.GetNumSeeds(); i++ {
		seedSet.Add(uint(rcQuery.GetSeed(i)))
	}
	for _, index := range matchingRCIndices {
		matchSet := m.index.GetSeedSet(index)
		if matchSet.CountIntersectionTo(seedSet, minRCMatches) < uint(minRCMatches) {
			continue
		}
		match := m.index.GetSeedSequence(index)
		seedMatches := match.Match(rcQuery, seedSet, matchSet, minRCMatches, k)
		if seedMatches != nil {
			for _, seedMatch := range seedMatches {
				start := match.GetOffset() + match.GetSeedOffset(seedMatch.MatchB[0], k)
				//qOffset/inset are swapped as they are based on the reverse-complement query, not the original query
				end := m.reference.Len() - match.GetInset() - match.GetSeedOffsetFromEnd(seedMatch.MatchB[len(seedMatch.MatchB)-1], k)
				if m.circular && start > m.reference.Len() {
					start -= m.reference.Len() //in the RHS of the circular join
				}
				qInset := rcQuery.GetSeedOffset(seedMatch.MatchA[0], k)
				qOffset := rcQuery.GetSeedOffsetFromEnd(seedMatch.MatchA[len(seedMatch.MatchA)-1], k)
				if qOffset+qInset > (rcQuery.Len()*2)/3 {
					continue
				}
				qInset += rcQuery.GetOffset()
				qOffset += rcQuery.GetInset()
				_, ids := seedMatch.GetBasesCovered(k)
				results = append(results, &Mapping{Start: start, End: end, QueryOffset: qOffset, QueryInset: qInset, RC: true, match: seedMatch, ids: ids})
				limit := (len(seedMatch.MatchA) * 4) / 5
				if limit > minRCMatches {
					minRCMatches = limit
				}
			}
		}
	}
	if len(results) > 1 {
		//sort and remove duplicate mappings
		sort.Sort(mappingsByPos(results))
		for i := len(results) - 1; i > 0; i-- {
			ra := results[i-1]
			rb := results[i]
			if ra.RC == rb.RC && rb.Start < ra.End {
				//keep the longer match
				if ra.End-ra.Start > rb.End-rb.Start {
					results[i] = results[len(results)-1]
					results = results[:len(results)-1]
				} else {
					results[i-1] = results[i]
					results[i] = results[len(results)-1]
					results = results[:len(results)-1]
				}
			}
		}
	}

	return results
}

func (m *mapper) MapWorker(queries <-chan sequence.Sequence, results chan<- []*Mapping, done chan<- bool) {
	for query := range queries {
		results <- m.Map(query)
	}
	done <- true
}
