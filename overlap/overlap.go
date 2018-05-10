package overlap

import (
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/seeds/alignment"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/util"
	"fmt"
)

//SeedQuery describes the input used to produce Overlaps from the current state of the Overlapper
type SeedQuery struct {
	ID                int //shared between forward and reverse-complement queries
	SequenceID        int
	Query             *seeds.SeedSequence
	AtStart           bool //adjacent to the beginning (if not, the end)
	ReverseComplement bool
}

var QueryEdges = 1
var QueryCentre = 2
var QueryAll = 4
var WeightEdges = 8
var WeightNone = 0 //default weighting

type Overlapper interface {
	PrepareQueries(int, int, []float64, <-chan sequence.Sequence, int) []*SeedQuery
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

func addWeighted(subseq sequence.Sequence, subseqsOut chan<- sequence.Sequence) {
	sideSize := 200
	if subseq.Len() > 400 {
		subseqsOut <- subseq.SubSequence(0, sideSize)
		subseqsOut <- subseq.SubSequence(subseq.Len()-sideSize, subseq.Len())
	} else {
		subseqsOut <- subseq
	}
}

func (lap *overlapper) getEdges(numSeeds, seedLimit int, weightSides bool, seqsIn <-chan sequence.Sequence, subseqsOut chan<- sequence.Sequence) []sequence.Sequence {
	cached := make([]sequence.Sequence, 0, 5000)
	for s := range seqsIn {
		if lap.index.Size() >= seedLimit {
			break
		}
		if s.Len() < lap.overlap*2 {
			if weightSides {
				addWeighted(s, subseqsOut)
			} else {
				subseqsOut <- s
			}
			cached = append(cached, s)
		} else {
			//fmt.Printf(">%s\n%s\n",s.GetName(),s.String())
			s1 := s.SubSequence(0, lap.overlap)
			s2 := s.SubSequence(s.Len()-lap.overlap, s.Len())
			if weightSides {
				addWeighted(s1, subseqsOut)
				addWeighted(s2, subseqsOut)
			} else {
				subseqsOut <- s1
				subseqsOut <- s2
			}
			cached = append(cached, s1)
			cached = append(cached, s2)
			//fmt.Printf(">%s_front\n%s\n",s.GetName(),s1.String())
			//fmt.Printf(">%s_back\n%s\n",s.GetName(),s2.String())
		}
	}
	//drain the remaining input, discarding the sequences
	for _ = range seqsIn {
	}
	return cached
}

func (lap *overlapper) getCentres(numSeeds, seedLimit int, weightSides bool, seqsIn <-chan sequence.Sequence, subseqsOut chan<- sequence.Sequence) []sequence.Sequence {
	cached := make([]sequence.Sequence, 0, 5000)
	for s := range seqsIn {
		if lap.index.Size() >= seedLimit {
			break
		}
		start := (s.Len() - lap.overlap) / 2
		if start < 0 {
			start = 0
		}
		end := start + lap.overlap
		if end >= s.Len() {
			end = s.Len() - 1
		}
		centre := s.SubSequence(start, end)
		if weightSides {
			addWeighted(centre, subseqsOut)
		} else {
			subseqsOut <- centre
		}
		cached = append(cached, centre)
	}
	//drain the remaining input, discarding the sequences
	for _ = range seqsIn {
	}
	return cached
}

func (lap *overlapper) getAll(numSeeds, seedLimit int, weightSides bool, seqsIn <-chan sequence.Sequence, subseqsOut chan<- sequence.Sequence) []sequence.Sequence {
	cached := make([]sequence.Sequence, 0, 5000)
	for s := range seqsIn {
		if lap.index.Size() >= seedLimit {
			break
		}
		fmt.Println("next sequence has length",s.Len())
		if s.Len() < lap.overlap*2 {
			if weightSides {
				addWeighted(s, subseqsOut)
			} else {
				subseqsOut <- s
			}
			cached = append(cached, s)
		} else {
			slices := s.Len() / lap.overlap
			for i := 0; i < slices; i++ {
				start := (i * s.Len())/slices
				end := ((i+1) * s.Len())/slices
				if i == slices-1 {
					end = s.Len()
				}
				sub := s.SubSequence(start,end)
				if weightSides {
					addWeighted(sub, subseqsOut)
				} else {
					subseqsOut <- sub
				}
				cached = append(cached, sub)
			}
		}
	}
	//drain the remaining input, discarding the sequences
	for _ = range seqsIn {
	}
	return cached
}

func (lap *overlapper) PrepareQueries(numSeeds int, seedLimit int, kmerValues []float64, seqs <-chan sequence.Sequence, queryType int) []*SeedQuery {
	completed := make(chan bool, lap.numWorkers)
	inputSeq := make(chan sequence.Sequence, lap.numWorkers*2)
	weightSides := (queryType & WeightEdges) != 0
	if weightSides {
		numSeeds /= 2
	}
	for i := 0; i < lap.numWorkers; i++ {
		go seeds.AddSeedsWorker(inputSeq, lap.index, numSeeds, kmerValues, completed)
	}
	var cached []sequence.Sequence
	if (queryType & QueryEdges) != 0 {
		cached = lap.getEdges(numSeeds, seedLimit, weightSides, seqs, inputSeq)
	} else if (queryType & QueryCentre) != 0 {
		cached = lap.getCentres(numSeeds, seedLimit, weightSides, seqs, inputSeq)
	} else { //query all default
		cached = lap.getAll(numSeeds, seedLimit, weightSides, seqs, inputSeq)
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
	for i, s := range cached {
		q := SeedQuery{queryID, s.GetID(), lap.index.NewSeedSequence(s), true, false}
		queries = append(queries, &q)
		rc := SeedQuery{queryID, q.SequenceID, q.Query.ReverseComplement(k,lap.index), true, true}
		queryID++
		queries = append(queries, &rc)
		cached[i] = nil
	}
	return queries
}

//AddSequences takes a channel of sequences and reads of a large set to add edges of for querying against.
func (lap *overlapper) AddSequences(seqs <-chan sequence.Sequence) {
	completed := make(chan bool, lap.numWorkers)
	chunksDone := make(chan bool, lap.numWorkers)
	seedSeq := make(chan *seeds.SeedSequence, lap.numWorkers*2)
	inputSeq := make(chan sequence.Sequence, lap.numWorkers*2)

	//start workers for getting seeds from sequences
	for i := 0; i < lap.numWorkers; i++ {
		go seeds.AddSequenceWorker(inputSeq, lap.index, seedSeq, completed)
	}
	//and our own set of workers to add them (split into chunks) as they complete
	for i := 0; i < lap.numWorkers; i++ {
		go lap.chunkWorker(seedSeq, chunksDone)
	}

	//then feed them in
	go func() {
		for s := range seqs {
			if s == nil {
				continue
			}
			inputSeq <- s
		}
		close(inputSeq)
		for i := 0; i < lap.numWorkers; i++ {
			<-completed
		}
		close(seedSeq)
	}()
	for i := 0; i < lap.numWorkers; i++ {
		<-chunksDone
	}
	lap.index.IndexSequences(lap.numWorkers)
}


func (lap *overlapper) chunkWorker(seedSeq <-chan *seeds.SeedSequence, done chan<- bool) {
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
			for {
				//count up seeds until 100 or past chunkSize bases
				seedCount := 0
				if prevSeedIndex >= s.GetNumSeeds()-150 { //we will add right up to the end
					if prevSeedIndex == 0 { //no need for a subsequence. Should have been caught as numChunks==1
						lap.index.AddSequence(s)
					} else {
						newFirstGap := s.GetNextSeedOffset(prevSeedIndex-1, k) - k
						lengthInBases += s.GetSeedOffsetFromEnd(prevSeedIndex, k) + k + newFirstGap
						lap.index.AddSequence(s.SubSequence(prevSeedIndex, s.GetNumSeeds()-1, lengthInBases, totalOffset-newFirstGap, 0))
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
				} else {
					//didn't get the minimum number of seeds, but we're up to chunkSize. Ignore this chunk entirely -- it contains no queries.
					prevSeedIndex += seedCount
					//move back overlap/2
					for seedCount = 0; lengthInBases < lap.overlap/2 && prevSeedIndex > 0; seedCount++ {
						prevSeedIndex--
						step := s.GetNextSeedOffset(prevSeedIndex, k)
						lengthInBases += step
						totalOffset -= step
					}
					lengthInBases = 0
				}
			}
		}
	}
	done <- true
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
	aligner := alignment.NewSeedAligner(lap.overlap/2)
	for q := range input {
		seedSet.Clear()
		for i := 0; i < q.Query.GetNumSeeds(); i++ {
			seedSet.Add(uint(q.Query.GetSeed(i)))
		}
		matches := lap.index.Matches(q.Query, lap.hitFraction)
		minMatches := int(lap.hitFraction*float64(q.Query.GetNumSeeds()) + 0.5)
		for _, match := range matches {
			matchSet := lap.index.GetSeedSet(match)
			if matchSet.CountIntersectionTo(seedSet, minMatches) < uint(minMatches) {
				continue
			}
			m := lap.index.GetSeedSequence(match)
			//TODO: if there are a vast number of matches (a multi-repeat) then move to a faster check here? Downstream won't be able to make much use of these anyway. Perhaps stick to strictly best match within this query, or first good match?
			//match is an index in the seed index, the sequence ID is external
			sMatches := aligner.PairwiseAlignments(q.Query, m, seedSet, matchSet, minMatches, k,false)
			if sMatches != nil {
				//choose the best one
				var best *seeds.SeedMatch
				bestCount := 0
				for _, sMatch := range sMatches {
					_, c := sMatch.GetBasesCovered(k)
					if c > bestCount {
						best = sMatch
					}
				}
				best.QueryID = q.ID
				best.ReverseComplementQuery = q.ReverseComplement
				output <- best
				//and some later matches can be pruned here too
				if len(best.MatchA)*2 > minMatches*3 {
					minMatches = (len(best.MatchA)*2)/3
				}
			}
		}
	}
	done <- true
}
