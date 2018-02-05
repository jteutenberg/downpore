package mapping

import (
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/util"
)

type Mapping struct {
	Start int
	End   int
	RC    bool
	match *seeds.SeedMatch
}

type Mapper interface {
	Map(sequence.Sequence) []*Mapping
}

type mapper struct {
	index    *seeds.SeedIndex
	edgeSize int //number of bases to use in queries
}

func NewMapper(reference sequence.Sequence, k uint, kmerValues []float64, seedRate int, numWorkers int) Mapper {
	m := mapper{index: seeds.NewSeedIndex(k), edgeSize: 1000}
	//walk the sequence adding 1 seed per seedRate bases
	m.index.AddSeeds(reference, seedRate, kmerValues)
	//index chunks of the reference
	chunks := make(chan sequence.Sequence, numWorkers*2)
	seedSeqs := make(chan *seeds.SeedSequence, numWorkers*2)
	done := make(chan bool, numWorkers)
	for i := 0; i < numWorkers; i++ {
		go seeds.AddSequenceWorker(chunks, m.index, seedSeqs, done)
	}
	ind := 0
	chunkSize := 10000
	overlapSize := 1000
	go func() {
		for j := 0; j < 10; j++ {
			start := j * chunkSize
			step := chunkSize*10 - overlapSize
			for i := start; i < reference.Len()-chunkSize; i += step {
				end := i + chunkSize
				if i >= reference.Len() {
					end = reference.Len()
				}
				//TODO: handle the short sequence at the end?
				//TODO: circular genomes
				chunks <- reference.SubSequence(i, end)
			}
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

func (m *mapper) Map(query sequence.Sequence) []*Mapping {
	var results []*Mapping
	//For short, return the front of the matches (ignoring substantially worse ones)
	if query.Len() <= m.edgeSize*2 {
		results = m.performMapping(query)
		//TODO: sort results from best to worst (by bases covered?)
	} else {
		//slice off the two ends of long sequences
		resultsA := m.performMapping(query.SubSequence(0, m.edgeSize))
		resultsB := m.performMapping(query.SubSequence(query.Len()-m.edgeSize, query.Len()))
		//For long, find reasonable pairings of start/end matches (if they exist)
		//pairs := make([]int, 0, len(resultsA)+len(resultsB))
		for _, ra := range resultsA {
			for _, rb := range resultsB {
				if ra.RC != rb.RC {
					continue //TODO: test for chimeric mapping here?
				}
				if ra.RC {
					ra, rb = rb, ra
				}
				distance := rb.Start - ra.End
				expectedDistance := ra.match.SeqA.Len() - m.edgeSize*2
				if expectedDistance < (distance*3)/2 && expectedDistance > (distance*2)/3 {
					//close enough? Longer distances really ought to be closer than 2/3rd, so tend towards 4/5?
					results = append(results, &Mapping{Start: ra.Start, End: rb.End, RC: ra.RC})
				}
			}
		}
	}
	return results
}

func (m *mapper) performMapping(query sequence.Sequence) []*Mapping {
	k := int(m.index.GetSeedLength())
	seedQuery := m.index.NewSeedSequence(query)
	rcQuery := seedQuery.ReverseComplement(k)
	matchingIndices := m.index.Matches(seedQuery, 0.3)
	matchingRCIndices := m.index.Matches(rcQuery, 0.3)
	minMatches := seedQuery.GetNumSeeds() / 3

	results := make([]*Mapping, 0, len(matchingIndices)+len(matchingRCIndices))

	seedSet := util.NewIntSet()
	//find the best hit
	for i := 0; i < seedQuery.GetNumSeeds(); i++ {
		seedSet.Add(uint(seedQuery.GetSeed(i)))
	}
	for _, index := range matchingIndices {
		match := m.index.GetSeedSequence(index)
		seedMatch := match.Match(seedQuery, seedSet, minMatches, k)
		if seedMatch != nil {
			start := int(index) * 2500
			//remove difference in offset between a and b
			start -= seedQuery.GetSeedOffset(seedMatch.MatchA[0], k) - match.GetSeedOffset(seedMatch.MatchB[0], k)
			end := start + query.Len() //TODO: wrong. Need to recalculate match end from last seed
			results = append(results, &Mapping{Start: start, End: end, RC: false, match: seedMatch})
		}
	}
	//or reverse-complement match
	seedSet.Clear()
	for i := 0; i < seedQuery.GetNumSeeds(); i++ {
		seedSet.Add(uint(rcQuery.GetSeed(i)))
	}
	for _, index := range matchingRCIndices {
		match := m.index.GetSeedSequence(index)
		seedMatch := match.Match(rcQuery, seedSet, minMatches, k)
		if seedMatch != nil {
			start := int(index) * 2500
			//remove difference in offset between a and b
			start -= rcQuery.GetSeedOffset(seedMatch.MatchA[0], k) - match.GetSeedOffset(seedMatch.MatchB[0], k)
			end := start + query.Len()
			results = append(results, &Mapping{Start: start, End: end, RC: true, match: seedMatch})
		}
	}

	return results
}
