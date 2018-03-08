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
	Query sequence.Sequence
	Start int
	End   int
	QueryOffset int
	QueryInset int
	RC    bool
	match *seeds.SeedMatch
	ids   int //absolute count of identity matches
}

type Mapper interface {
	Map(sequence.Sequence) []*Mapping
	MapWorker(<-chan sequence.Sequence, chan<- []*Mapping, chan<- bool)
	AsString(*Mapping) string
}

type mapper struct {
	index    *seeds.SeedIndex
	reference sequence.Sequence
	edgeSize int //number of bases to use in queries
	circular bool
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
func (m mappingsByPos) Swap(i,j int) {
	m[i],m[j] = m[j],m[i]
}
func (m mappingsByQPos) Len() int {
	return len(m)
}
func (m mappingsByQPos) Less(i, j int) bool {
	return m[i].QueryOffset < m[j].QueryOffset
}
func (m mappingsByQPos) Swap(i,j int) {
	m[i],m[j] = m[j],m[i]
}
func (m mappingsByIds) Len() int {
	return len(m)
}
func (m mappingsByIds) Less(i, j int) bool {
	return m[i].ids > m[j].ids
}
func (m mappingsByIds) Swap(i,j int) {
	m[i],m[j] = m[j],m[i]
}

func NewMapper(reference sequence.Sequence, circular bool, k uint, kmerValues []float64, seedRate int, edgeSize int, numWorkers int) Mapper {
	m := mapper{index: seeds.NewSeedIndex(k), reference:reference, edgeSize: edgeSize, circular:circular}
	//walk the sequence adding 1 seed per seedRate bases
	m.index.AddSingleSeeds(reference, seedRate, kmerValues)
	log.Println("Index has",m.index.Size(),"seeds")
	//index chunks of the reference
	chunks := make(chan sequence.Sequence, numWorkers*2)
	seedSeqs := make(chan *seeds.SeedSequence, numWorkers*2)
	done := make(chan bool, numWorkers)
	for i := 0; i < numWorkers; i++ {
		go seeds.AddSequenceWorker(chunks, m.index, seedSeqs, done)
	}
	ind := 0
	chunkSize := 10000
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
			chunks <- reference.SubSequence(reference.Len()-edgeSize, reference.Len()).Append(0,reference.SubSequence(0,edgeSize),nil)
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
	return &m
}

//AsString prepares a PAF line for this mapping
func (m *mapper) AsString(mapping *Mapping) string {
	rc := "+"
	if mapping.RC {
		rc = "-"
	}
	return fmt.Sprintf("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t255",mapping.Query.GetName(),mapping.Query.Len(),mapping.QueryOffset,mapping.Query.Len()-mapping.QueryInset,rc,m.reference.GetName(),m.reference.Len(),mapping.Start,mapping.End,mapping.ids,mapping.End-mapping.Start)
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
		distance = right.Start - left.End// - qInset + inset - qOffset + offset //TODO: <- problem: inset/offset change over time
	} else {
		distance = left.Start - right.End// - qInset + inset - qOffset + offset
	}
	if m.circular && distance < -50 {
		distance += m.reference.Len()
	}
	return (distance < 50 && expectedDistance < 50 && distance > -50) || (expectedDistance < (distance*3)/2 && expectedDistance > (distance*2)/3)
}

//maps the (potentially noisy) ends
//if any map end-to-end then these are merged and all others discarded
func (m *mapper) mapEnds(query sequence.Sequence) (openA, openB, matching []*Mapping) {
	openA = m.performMapping(query.SubSequence(0, m.edgeSize))
	openB = m.performMapping(query.SubSequence(query.Len()-m.edgeSize, query.Len()))
	openA = removeDominated(openA,openA,query.Len())
	openB = removeDominated(openB,openB,query.Len())
	updateQuery(openA,query)
	updateQuery(openB,query)
	return m.matchPairs(openA, openB)
}

func (m *mapper) matchPairs(openA, openB []*Mapping) (remainingA, remainingB, matched []*Mapping) {
	for i := len(openA)-1; i >= 0; i-- {
		ra := openA[i]
		for j := len(openB)-1; j >= 0; j-- {
			rb := openB[j]
			if m.isConsistent(ra, rb) {
				qOffset := ra.QueryOffset
				qInset := rb.QueryInset
				if ra.RC {
					ra, rb = rb, ra
				}
				combined := &Mapping{Start: ra.Start, End: rb.End, Query:ra.Query, QueryOffset:qOffset, QueryInset:qInset, RC: ra.RC, ids:ra.ids+rb.ids}
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

func (m *mapper) GetRepeats(start, end int) (starts, ends []int) {
	query := m.reference.SubSequence(start,end)
	hits := m.performMapping(query)
	hits = removeDominated(hits,hits,query.Len())
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
		newA = removeDominated(newA,newA,query.Len())
		updateQuery(newA,query)
		openA, newA, extended = m.matchPairs(openA, newA)
		if extended != nil {
			openA = append(newA, extended...)
		} else {
			openA = append(openA,newA...)
		}
		//match with B. If nothing works then this is a nonsense sequence
		newA, newB, matched = m.matchPairs(openA, openB)
		if matched == nil {
			return newA[:0],newB[:0],extended
		}
		return newA[:0],newB[:0],matched
	}

	// Long sequences

	// 1. Test one in from the edges, see if anything matches up.
	newA = m.performMapping(query.SubSequence(m.edgeSize, m.edgeSize*2))
	newA = removeDominated(newA,newA,query.Len())
	updateQuery(newA,query)
	debug := false
	if debug {
		for i, a := range openA {
			log.Println(" -a- ",i,": ",a.Start,a.End," :: ",a.QueryOffset,a.Query.Len()-a.QueryInset)
		}
		for i, a := range newA {
			log.Println(" -+a- ",i,": ",a.Start,a.End)
		}
	}
	openA, newA, extended = m.matchPairs(openA, newA)
	if extended == nil {
		//log.Println("No matches in A!")
		openA = newA //drop the edges, but keep looking for matches
	} else {
		if debug {log.Println(len(extended),"extensions in A:")
			for i, a := range extended {
				log.Println(" -xa- ",i,": ",a.Start,a.End)
			}
		}
		openA = append(newA,extended...)
	}
	//do the same from the other end
	newB = m.performMapping(query.SubSequence(query.Len()-m.edgeSize*2, query.Len()-m.edgeSize))
	newB = removeDominated(newB,newB,query.Len())
	updateQuery(newB,query)
	if debug {
		for i, a := range openB {
			log.Println(" -b- ",i,": ",a.Start,a.End)
		}
		for i, a := range newB {
			log.Println(" -+b- ",i,": ",a.Start,a.End)
		}
	}
	openB, newB, extended = m.matchPairs(newB, openB)
	if extended == nil {
		//log.Println("No matches in B!")
		openB = newB
	} else {
		if debug {
			log.Println(len(extended),"extensions in B:")
			for i, a := range extended {
				log.Println(" -xb- ",i,": ",a.Start,a.End)
			}
		}
		openB = append(newB, extended...)
	}

	// see if anything matches full length now
	newA, newB, matched = m.matchPairs(openA, openB)

	// 2. Second (last!) chance: take two steps in from each edge and try again
	if matched == nil {
		//do a second round, just in case we got unlucky with read errors or reference seeds
		if query.Len() > m.edgeSize*6 {
			openA = m.performMapping(query.SubSequence(m.edgeSize*2, m.edgeSize*3))
			openA = removeDominated(openA,openA,query.Len())
			updateQuery(openA, query)
			if debug {
				for i, a := range openA {
					log.Println(" -++a- ",i,": ",a.Start,a.End)
				}
			}
			openA, newA, extended = m.matchPairs(newA, openA)
			if extended != nil {
				if debug {
					log.Println(len(extended),"extensions in A:")
					for i, a := range extended {
						log.Println(" -xxa- ",i,": ",a.Start,a.End)
					}
				}
				openA = append(openA,extended...)
			} else if debug {
				log.Println("No extensions in A")
			}
			openA = append(openA,newA...) //keep everything for now. One last chance to match with B
			openB = m.performMapping(query.SubSequence(query.Len()-m.edgeSize*3, query.Len()-m.edgeSize*2))
			openB = removeDominated(openB,openB,query.Len())
			updateQuery(openB, query)
			if debug {
				for i, a := range openB {
					log.Println(" -++b- ",i,": ",a.Start,a.End)
				}
			}
			openB, newB, extended = m.matchPairs(openB,newB)
			if extended != nil {
				if debug {
					log.Println(len(extended),"extensions in B:")
					for i, a := range extended {
						log.Println(" -xxb- ",i,": ",a.Start,a.End)
					}
				}
				openB = append(openB, extended...)
			} else if debug {
				log.Println("No extensions in B")
			}
			openB = append(openB, newB...)
			newA, newB, matched = m.matchPairs(openA,openB)
		}
	}
	//log.Println("Final sizes:",len(newA),len(newB),len(matched))
	return newA, newB, matched
}

//removeDominated modifies the ordering of open, returning a new slice containing only those mappings not dominated by a member of extended.
//The extended slice must already be ordered by query position
func removeDominated(open, extended []*Mapping, queryLen int) ([]*Mapping) {
	if len(open) == 0 || len(extended) == 0 {
		return open
	}
	sort.Sort(mappingsByQPos(open))
	j := 0 //first relevant index in extended
	//90% contained, 20% less identity match
	toRemove := make([]bool, len(open))
	for i,next := range open {
		for j < len(extended) && queryLen - extended[j].QueryInset < next.QueryOffset {
			j++
		}
		if j == len(extended) {
			return open
		}
		//test all relevant in extended
		dominated := false
		for k := j; !dominated && k < len(extended) && extended[k].QueryOffset < queryLen - next.QueryInset; k++ {
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
				dominated = ((end-start)*10 > (queryLen-next.QueryOffset-next.QueryInset)*9 )
			}
		}
		toRemove[i] = dominated
	}
	last := len(open)-1
	for i := last; i >= 0; i-- {
		if toRemove[i] {
			open[i] = open[last]
			last--
		}
	}
	return open[:last+1]
}

//TODO: a function to search for chimeric boundary given no end-to-end matches, but some consistent on each side

//makes some additional mappings to ensure the query (already mapped up to left/right) is a full mapping
func (m *mapper) validateMapping(query sequence.Sequence, mapping *Mapping, mappedLeft, mappedRight int) bool {
	return true
}

func (m *mapper) Map(query sequence.Sequence) []*Mapping {
	var results []*Mapping
	//For short, return the front of the matches (ignoring substantially worse ones)
	if query.Len() <= m.edgeSize*2 {
		results = m.performMapping(query)
		results = removeDominated(results,results,query.Len())
		updateQuery(results,query)
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
			openA, openB, matched = m.mapNext(query,openA, openB)
			if matched != nil {
				//the end-to-end mappings are still the only mappings of interest
				results = matched
			} else {
				//Otherwise: remove all unpaired ends
				size := query.Len()-m.edgeSize
				for i := len(openA)-1; i >= 0; i-- {
					if openA[i].QueryInset >= size {
						openA[i] = openA[len(openA)-1]
						openA = openA[:len(openA)-1]
					}
				}
				for i := len(openB)-1; i >= 0; i-- {
					if openB[i].QueryOffset >= size {
						openB[i] = openB[len(openB)-1]
						openB = openB[:len(openB)-1]
					}
				}
				results = append(openA, openB...)
			}
		}
		//TODO: sort results from best to worst (by bases covered?)
	}
	/*if len(results) > 0 {
		log.Println("--")
		for _, r := range results {
			log.Println("q:",r.Start,r.End,r.RC)
			if !r.RC && r.Start < r.End {
				ss, es := m.GetRepeats(r.Start, r.End)
				for i := 0; i < len(ss); i++ {
					log.Println("rep:",ss[i],es[i])
				}
			}
		}
	}*/
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
	seedSet := util.NewIntSetCapacity(maxSeed+1)
	for i := 0; i < seedQuery.GetNumSeeds(); i++ {
		seedSet.Add(uint(seedQuery.GetSeed(i)))
	}
	for _, index := range matchingIndices {
		//1. Does this have enough matching seeds (compared to best so far)?
		matchSet := m.index.GetSeedSet(index)
		/*if index < 10 {
			log.Println(matchSet.CountIntersection(seedSet),"at index",index,"out of",matchSet.CountMembers(),"/",seedSet.CountMembers(),"and min is",minMatches)
		}*/
		if matchSet.CountIntersectionTo(seedSet,minMatches) < uint(minMatches) {
			continue
		}
		//2. Match based on shared seeds
		match := m.index.GetSeedSequence(index)
		seedMatches := match.Match(seedQuery, seedSet, matchSet, minMatches, k)
		if seedMatches != nil {
			for _, seedMatch := range seedMatches {
				/*if index < 10 {
					log.Println("Forward match on offset",seedQuery.GetOffset(),seedQuery.GetInset(),"to",match.GetOffset())
					log.Println(seedMatch.LongString(m.index))
				}*/
				start := match.GetOffset() + match.GetSeedOffset(seedMatch.MatchB[0],k)
				end := m.reference.Len()-match.GetInset()-match.GetSeedOffsetFromEnd(seedMatch.MatchB[len(seedMatch.MatchB)-1],k)
				qOffset := seedQuery.GetSeedOffset(seedMatch.MatchA[0], k)
				qInset := seedQuery.GetSeedOffsetFromEnd(seedMatch.MatchA[len(seedMatch.MatchA)-1],k)
				if qOffset+qInset > (seedQuery.Len()*2)/3 {
					continue
				}
				qOffset += seedQuery.GetOffset()
				qInset += seedQuery.GetInset()
				_, ids := seedMatch.GetBasesCovered(k) //bases matched in the reference
				results = append(results, &Mapping{Start: start, End: end, QueryOffset:qOffset, QueryInset:qInset,RC: false, match: seedMatch, ids:ids })
				/*if index < 10 {
					log.Println(start,end,"from",qOffset,qInset)
				}*/
				limit := (len(seedMatch.MatchA)*4)/5
				if limit > minMatches {
					minMatches = limit
				}
				if limit > minRCMatches {
					minRCMatches = limit
				}
			} /*else if index < 10 {
				log.Println("No seed match though.")
			}*/
		}
	}
	//log.Println(query.GetName(),"Num seeds:",seedQuery.GetNumSeeds(),rcQuery.GetNumSeeds(),"unique:",seedSet.CountMembers())
	//log.Println("matched",len(results),"/",len(matchingIndices))
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
				//log.Println("Reverse match on offset",rcQuery.GetOffset(),rcQuery.GetInset(),"to",match.GetOffset())
				//log.Println(seedMatch.LongString(m.index))
				start := match.GetOffset() + match.GetSeedOffset(seedMatch.MatchB[0],k)
				//qOffset/inset are swapped as they are based on the reverse-complement query, not the original query
				end := m.reference.Len()-match.GetInset()-match.GetSeedOffsetFromEnd(seedMatch.MatchB[len(seedMatch.MatchB)-1],k)
				qInset := rcQuery.GetSeedOffset(seedMatch.MatchA[0], k)
				qOffset := rcQuery.GetSeedOffsetFromEnd(seedMatch.MatchA[len(seedMatch.MatchA)-1],k)
				if qOffset+qInset > (rcQuery.Len()*2)/3 {
					continue
				}
				qInset += rcQuery.GetOffset()
				qOffset += rcQuery.GetInset()
				_, ids := seedMatch.GetBasesCovered(k)
				results = append(results, &Mapping{Start: start, End: end, QueryOffset: qOffset, QueryInset:qInset, RC: true, match: seedMatch, ids:ids})
				limit := (len(seedMatch.MatchA)*4)/5
				if limit > minRCMatches {
					minRCMatches = limit
				}
			}
		}
	}
	if len(results) > 1 {
		//sort and remove duplicate mappings
		sort.Sort(mappingsByPos(results))
		for i := len(results)-1; i > 0; i-- {
			ra := results[i-1]
			rb := results[i]
			if ra.RC == rb.RC && rb.Start < ra.End {
				//keep the longer match
				if ra.End - ra.Start > rb.End - rb.Start {
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
