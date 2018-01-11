package seeds

import (
	"fmt"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/util"
	"math"
	"sync"
)

//SeedIndex is a set of reads in gapped-seed format ready for querying
type SeedIndex struct {
	seedSize     uint
	seeds        []bool
	sequences    []*SeedSequence //list of all sequences
	sequenceSets []*util.IntSet  //list of seed -> set of sequences (indices)
	size         uint
	lock         *sync.Mutex
}

func NewSeedIndex(k uint) *SeedIndex {
	size := 1
	for j := k; j > 0; j-- {
		size *= 4
	}
	var lock sync.Mutex
	sg := SeedIndex{seeds: make([]bool, size, size), seedSize: k, sequences: make([]*SeedSequence, 0, 100000), sequenceSets: make([]*util.IntSet, size, size), size: 0, lock: &lock}
	return &sg
}

//NewSeedSequence re-uses any seeds in the index, adding additional seeds as required to
//reach at least minSeeds for the sequence. Note that while seeds are added to the index, the
//resulting sequence itself is not added.
func (g *SeedIndex) NewSeedSequence(seq sequence.Sequence, minSeeds int, ranks []float64) *SeedSequence {
	//first-pass, get the sequence as kmers
	kmers := seq.Kmers(int(g.seedSize))
	//count the number existing in the index and generate the gapped-seed sequence using those
	count := 0
	segments := make([]int, 0, 20)
	prev := 0
	for i, seed := range kmers {
		if g.seeds[seed] {
			count++
			segments = append(segments, i-prev)
			segments = append(segments, int(seed))
			prev = i + int(g.seedSize)
		}
	}
	segments = append(segments, len(kmers)-prev+int(g.seedSize)-1) //number of bases remaining
	if count < minSeeds {
		topN := make([]uint, minSeeds-count, minSeeds-count)
		topNValues := make([]float64, minSeeds-count, minSeeds-count)
		for i := 0; i < len(topNValues); i++ {
			topNValues[i] = math.MaxFloat64
		}
		//add additional seeds, not overlapping existing seeds
		prevIndex := 0
		for i := 0; i < len(segments); i += 2 {
			//earliest and latest k-mer indices that can be seeds
			earliest := prevIndex
			latest := prevIndex + segments[i] - int(g.seedSize) //we'll do batches of k k-mers
			if i == len(segments)-1 {
				latest -= int(g.seedSize) //turn bases into k-mers (a shorter sequence)
			}
			for earliest < latest {
				//pick one from the first k k-mers, then step past it
				bestIndex := -1 //index offset from prevIndex
				bestValue := math.MaxFloat64
				//run over k k-mers
				for j := 0; j < int(g.seedSize); j++ {
					seed := kmers[earliest+j]
					value := ranks[seed]
					if value < bestValue {
						bestValue = value
						bestIndex = earliest + j
					}
				}
				//step past the one we chose (the best)
				earliest = bestIndex + int(g.seedSize)
				//and also add it to our collection
				n := 0
				for ; n < len(topNValues) && topNValues[n] > bestValue; n++ {
					if n > 0 {
						topNValues[n-1] = topNValues[n]
						topN[n-1] = topN[n]
					}
				}
				if n > 0 {
					topNValues[n-1] = bestValue
					topN[n-1] = uint(kmers[bestIndex])
				}
			}
			prevIndex += segments[i] + int(g.seedSize)
		}
		// store the top seeds (synchronise here)
		g.lock.Lock()
		for _, seed := range topN {
			g.seeds[seed] = true
			g.size++
			seed = ReverseComplement(seed, g.seedSize)
			if !g.seeds[seed] {
				g.seeds[seed] = true
				g.size++
			}
		}
		g.lock.Unlock()

		//and regenerate the gapped-seed sequence
		prev = 0
		segments = segments[:0]
		for i, seed := range kmers {
			if g.seeds[seed] {
				segments = append(segments, i-prev)
				segments = append(segments, int(seed))
				prev = i + int(g.seedSize)
			}
		}
		segments = append(segments, len(kmers)-prev+int(g.seedSize))
	}
	name := seq.GetName()
	seedSeq := SeedSequence{segments: segments, length: seq.Len(), id: seq.GetID(), name: &name, offset: seq.GetOffset(), inset: seq.GetInset(), rc: false}
	return &seedSeq
}

//NewAllSeedSequence creates a new SeedSequence adding every seed to the index
//Unlike NewSeedSequence this will not add the reversecomplent seeds to the index
func (g *SeedIndex) NewAllSeedSequence(seq sequence.Sequence) *SeedSequence {
	kmers := seq.Kmers(int(g.seedSize))
	segments := make([]int, 0, seq.Len()*2)
	prev := 0
	g.lock.Lock()
	for i, seed := range kmers {
		if !g.seeds[seed] {
			g.seeds[seed] = true
			g.size++
		}
		segments = append(segments, i-prev)
		segments = append(segments, int(seed))
		prev = i + int(g.seedSize)
	}
	g.lock.Unlock()
	segments = append(segments, 0)
	name := seq.GetName()
	seedSeq := SeedSequence{segments: segments, length: seq.Len(), id: seq.GetID(), name: &name, offset: seq.GetOffset(), inset: seq.GetInset(), rc: false}
	return &seedSeq
}

func (g *SeedIndex) Size() uint {
	return g.size
}

func (g *SeedIndex) GetSeedSequence(index uint) *SeedSequence {
	return g.sequences[index]
}

func (g *SeedIndex) GetSeedLength() uint {
	return g.seedSize
}

func (g *SeedIndex) GetNumSequences() uint {
	return uint(len(g.sequences))
}

func (g *SeedIndex) Contains(seed uint) bool {
	return g.seeds[seed]
}

func (g *SeedIndex) AddSequence(seq *SeedSequence) {
	index := uint(len(g.sequences))
	g.lock.Lock()
	g.sequences = append(g.sequences, seq)
	for i, seed := range seq.segments {
		if (i & 1) == 1 {
			if g.sequenceSets[seed] == nil {
				g.sequenceSets[seed] = util.NewIntSet()
				if !g.seeds[seed] {
					g.seeds[seed] = true
					g.size++
				}
			}
			g.sequenceSets[seed].Add(index)
		}
	}
	g.lock.Unlock()
}

func (g *SeedIndex) RemoveSequences() {
	g.sequences = g.sequences[:0]
	for _, seqs := range g.sequenceSets {
		if seqs != nil {
			seqs.Clear()
		}
	}
}

func (g *SeedIndex) Destroy() {
	g.sequences = nil
	g.sequenceSets = nil
}

func (g *SeedIndex) PrintSeeds() {
	for i, seqs := range g.sequenceSets {
		if seqs != nil {
			fmt.Println("Seed", i, "in", seqs.CountMembers())
		}
	}
}

//Matches finds all sequences that contain hitFraction of all seeds in the query
func (g *SeedIndex) Matches(query *SeedSequence, hitFraction float64) []uint {
	allSeedSets := make([]*util.IntSet, 0, len(query.segments)/2)
	prevSeed := -1 //remove some redundant seeds
	//maxSeqs := uint(len(g.sequences)/2)
	maxSeqs := uint(len(g.sequences))
	for i := 1; i < len(query.segments); i += 2 {
		seed := query.segments[i]
		adj := g.sequenceSets[seed]
		if seed != prevSeed && adj != nil && adj.Size() < maxSeqs {
			allSeedSets = append(allSeedSets, adj)
			prevSeed = seed
		}
	}
	if len(allSeedSets) < 5 { //not many usable seeds in the query!
		return make([]uint, 0, 0)
	}
	minCount := int(hitFraction*float64(len(allSeedSets)) + 0.5)
	return util.GetSharedIDs(allSeedSets, minCount)
}

//AddSequence processes all sequences using the provided index but not adding
//any new seeds, outputting the equivalent SeedSequence
func AddSequenceWorker(input <-chan sequence.Sequence, index *SeedIndex, output chan<- *SeedSequence, done chan<- bool) {
	for seq := range input {
		output <- index.NewSeedSequence(seq, 0, nil)
	}
	done <- true
}

//AddSeedsWorker processes all sequences, adding seeds to the provided index
func AddSeedsWorker(input <-chan sequence.Sequence, index *SeedIndex, numSeeds int, kmerValues []float64, done chan<- bool) {
	for seq := range input {
		index.NewSeedSequence(seq, numSeeds, kmerValues)
	}
	done <- true
}
