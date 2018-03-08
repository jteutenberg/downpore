package overlap

import (
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/util"
)

//SeedQuery describes the input used to produce Overlaps from the current state of the Overlapper
type SeedQuery struct {
	ID                int //shared between forward and reverse-complement queries
	SequenceID        int
	Query             *seeds.SeedSequence
	AtStart           bool //adjacent to the beginning (if not, the end)
	ReverseComplement bool
}

type Overlapper interface {
	PrepareQueries(int, int, []float64, <-chan sequence.Sequence, bool) []*SeedQuery
	AddSequences(<-chan sequence.Sequence)
	FindOverlaps([]*SeedQuery) <-chan *seeds.SeedMatch
	SetOverlapSize(int)
}

type overlapper struct {
	index       *seeds.SeedIndex
	chunkSize   uint //slice sequences into this size chunks for querying against
	numWorkers  int
	overlap     int
	hitFraction float64
	minSeeds    int
}

func NewOverlapper(index *seeds.SeedIndex, chunkSize uint, numWorkers int, overlap int, minSeeds int, hitFraction float64) Overlapper {
	ov := overlapper{index, chunkSize, numWorkers, overlap, hitFraction, minSeeds}
	return &ov
}

//TODO: make this a factory-style method (currently only uses "centre")
func (lap *overlapper) PrepareQueries(numSeeds int, seedLimit int, kmerValues []float64, seqs <-chan sequence.Sequence, centre bool) []*SeedQuery {
	completed := make(chan bool, lap.numWorkers)
	inputSeq := make(chan sequence.Sequence, lap.numWorkers*2)
	weightSides := false
	if weightSides {
		numSeeds /= 2
	}
	for i := 0; i < lap.numWorkers; i++ {
		go seeds.AddSeedsWorker(inputSeq, lap.index, numSeeds, kmerValues, completed)
	}
	//TODO: fill the cache once, according to parameters
	cached := make([]sequence.Sequence, 0, 5000)
	for s := range seqs {
		if lap.index.Size() >= seedLimit {
			break
		}
		if centre {
			start := (s.Len() - lap.overlap) / 2
			if start < 0 {
				start = 0
			}
			end := start + lap.overlap
			if end >= s.Len() {
				end = s.Len() - 1
			}
			if weightSides {
				edge := 200
				if edge >= end {
					edge = end
				}
				inputSeq <- s.SubSequence(start, start+edge)
				inputSeq <- s.SubSequence(end-edge, end)
			} else {
				inputSeq <- s.SubSequence(start, end)
			}
			cached = append(cached, s.SubSequence(start, end))
		} else { //use edges of the sequence instead
			//enqueue this sequence for processing into a query
			if weightSides {
				inputSeq <- s.SubSequence(0, 200)
				inputSeq <- s.SubSequence(s.Len()-200, s.Len())
				if s.Len() >= lap.overlap*2 {
					inputSeq <- s.SubSequence(lap.overlap-200, lap.overlap)
					inputSeq <- s.SubSequence(s.Len()-lap.overlap, s.Len()-lap.overlap+200)
					cached = append(cached, s.SubSequence(0, lap.overlap))
					cached = append(cached, s.SubSequence(s.Len()-lap.overlap, s.Len()))
				} else {
					cached = append(cached, s)
				}
			} else {
				if s.Len() < lap.overlap*2 {
					inputSeq <- s
					cached = append(cached, s)
					//fmt.Printf(">%s\n%s\n",s.GetName(),s.String())
				} else {
					s1 := s.SubSequence(0, lap.overlap)
					s2 := s.SubSequence(s.Len()-lap.overlap, s.Len())
					inputSeq <- s1
					inputSeq <- s2
					cached = append(cached, s1)
					cached = append(cached, s2)
					//fmt.Printf(">%s_front\n%s\n",s.GetName(),s1.String())
					//fmt.Printf(">%s_back\n%s\n",s.GetName(),s2.String())
				}
			}
		}
	}
	//drain the remaining input, discarding the sequences
	for _ = range seqs {
	}
	//finish up the workers
	close(inputSeq)
	for i := 0; i < lap.numWorkers; i++ {
		<-completed
	}

	//recalculate all queries using the full set of seeds
	queries := make([]*SeedQuery, 0, len(cached))
	seedSeq := make(chan *seeds.SeedSequence, lap.numWorkers*2)
	inputSeq = make(chan sequence.Sequence, lap.numWorkers*2)
	for i := 0; i < lap.numWorkers; i++ {
		go seeds.AddSequenceWorker(inputSeq, lap.index, seedSeq, completed)
	}
	//linear here: grab the short queries using the full set of seeds
	queryID := 0
	k := int(lap.index.GetSeedLength())
	for _, s := range cached {
		q := SeedQuery{queryID, s.GetID(), lap.index.NewSeedSequence(s), true, false}
		queries = append(queries, &q)
		rc := SeedQuery{queryID, q.SequenceID, q.Query.ReverseComplement(k), true, true}
		queryID++
		queries = append(queries, &rc)
	}
	return queries
}

//AddSequences takes a channel of sequences and reads of a large set to add edges of for querying against.
func (lap *overlapper) AddSequences(seqs <-chan sequence.Sequence) {
	completed := make(chan bool, lap.numWorkers)
	seedSeq := make(chan *seeds.SeedSequence, lap.numWorkers*2)
	inputSeq := make(chan sequence.Sequence, lap.numWorkers*2)

	//start workers for getting seeds from sequences
	for i := 0; i < lap.numWorkers; i++ {
		go seeds.AddSequenceWorker(inputSeq, lap.index, seedSeq, completed)
	}
	count := lap.index.GetNumSequences()
	go func() {
		for s := range seqs {
			if s == nil {
				continue
			}
			inputSeq <- s
			numChunks := s.Len()/int(lap.chunkSize) + 1
			count += uint(numChunks) //later we will chop to approximately this number of chunks
		}
		close(inputSeq)
		for i := 0; i < lap.numWorkers; i++ {
			<-completed
		}
		close(seedSeq)
	}()
	//and add them as they get done
	k := int(lap.index.GetSeedLength())
	for s := range seedSeq {
		//slice into chunkSize pieces
		//TODO: make this chop part of SeedSequence itself?
		numChunks := s.Len()/int(lap.chunkSize) + 1
		if numChunks == 1 || s.GetNumSeeds() < lap.minSeeds*3 {
			if s.GetNumSeeds() >= lap.minSeeds {
				lap.index.AddSequence(s)
			}
		} else {
			prevSeedIndex := 0                   //chop by 100 seeds
			totalOffset := s.GetSeedOffset(0, k) //offset to the first seed
			lengthInBases := 0                   //initial gap will be added to the length later
			count := 0
			for {
				//count up seeds until 100 or past chunkSize bases
				seedCount := 0
				if prevSeedIndex >= s.GetNumSeeds()-150 { //we will add right up to the end
					if prevSeedIndex == 0 { //no need for a subsequence. Should have been caught as numChunks==1
						lap.index.AddSequence(s)
					} else {
						newFirstGap := s.GetNextSeedOffset(prevSeedIndex-1, k) - k
						lengthInBases += s.GetSeedOffsetFromEnd(prevSeedIndex, k) + k + newFirstGap
						//fmt.Println("total offset:",totalOffset,"less gap",newFirstGap,"plus length",lengthInBases,"should equal total length",s.GetLength(),". Note: actual seed offset is",s.GetSeedOffset(prevSeedIndex,k),"and offset from end is",s.GetSeedOffsetFromEnd(prevSeedIndex,k))
						lap.index.AddSequence(s.SubSequence(prevSeedIndex, s.GetNumSeeds()-1, lengthInBases, totalOffset-newFirstGap, 0)) //s.GetLength()-totalOffset+newFirstGap-lengthInBases))
					}
					break
				}
				for ; lengthInBases < int(lap.chunkSize) && seedCount < 100 && prevSeedIndex+seedCount < s.GetNumSeeds(); seedCount++ {
					lengthInBases += s.GetNextSeedOffset(prevSeedIndex+seedCount, k)
				}
				//perform the chop
				if seedCount >= lap.minSeeds {
					newFirstGap := s.GetNextSeedOffset(prevSeedIndex-1, k) - k
					lengthInBases += newFirstGap //add the gap, and the bases in the first seed
					lap.index.AddSequence(s.SubSequence(prevSeedIndex, prevSeedIndex+seedCount-1, lengthInBases, totalOffset-newFirstGap, s.GetLength()-totalOffset-lengthInBases+newFirstGap))
					totalOffset += lengthInBases - newFirstGap //this moves length in bases up to the following seed
					lengthInBases = 0
					prevSeedIndex += seedCount
					if prevSeedIndex >= s.GetNumSeeds() {
						break
					}
					//move back 5 seeds or overlap/2, whichever is larger
					for seedCount = 0; seedCount < 5 && lengthInBases < lap.overlap/2 && prevSeedIndex > 0; seedCount++ {
						prevSeedIndex--
						step := s.GetNextSeedOffset(prevSeedIndex, k)
						lengthInBases += step
						totalOffset -= step
					}
					lengthInBases = 0
					count++
				}
			}
		}
	}
}

func (lap *overlapper) FindOverlaps(queries []*SeedQuery) <-chan *seeds.SeedMatch {
	//start workers
	input := make(chan *SeedQuery, lap.numWorkers*2)
	done := make(chan bool, lap.numWorkers)
	output := make(chan *seeds.SeedMatch, lap.numWorkers*2)
	for i := 0; i < lap.numWorkers; i++ {
		go lap.matchWorker(input, output, done)
	}
	//feed in the queries
	go func() {
		for _, q := range queries {
			input <- q
		}
		close(input)
		for i := 0; i < lap.numWorkers; i++ {
			<-done
		}
		close(output)
	}()
	return output
}

func (lap *overlapper) SetOverlapSize(size int) {
	lap.overlap = size
}

func (lap *overlapper) matchWorker(input <-chan *SeedQuery, output chan<- *seeds.SeedMatch, done chan<- bool) {
	k := int(lap.index.GetSeedLength())
	seedSet := util.NewIntSet()
	for q := range input {
		seedSet.Clear()
		for i := 0; i < q.Query.GetNumSeeds(); i++ {
			seedSet.Add(uint(q.Query.GetSeed(i)))
		}
		matches := lap.index.Matches(q.Query, lap.hitFraction)
		minMatches := int(lap.hitFraction*float64(q.Query.GetNumSeeds()) + 0.5)
		for _, match := range matches {
			m := lap.index.GetSeedSequence(match)
			//match is an index in the seed index, the sequence ID is external
			sMatches := m.Match(q.Query, seedSet, lap.index.GetSeedSet(match), minMatches, k)
			if sMatches != nil {
				for _, sMatch := range sMatches {
					sMatch.QueryID = q.ID
					sMatch.ReverseComplementQuery = q.ReverseComplement
					output <- sMatch
				}
			}
		}
	}
	done <- true
}
