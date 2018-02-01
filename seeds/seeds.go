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

func (g *SeedIndex) NewSeedSequence(seq sequence.Sequence) *SeedSequence {
	k := int(g.seedSize)
	//for speed, we duplicate Kmers() for a single pass
	var mask int
	for i := 0; i < k; i++ {
		mask = (mask << 2) | 3
	}
	count := seq.CountKmers(k, mask, g.seeds)
	segments := make([]int, count*2+1)
	seq.WriteSegments(segments, k, mask, g.seeds)
	name := seq.GetName()
	seedSeq := SeedSequence{segments: segments, length: seq.Len(), id: seq.GetID(), name: &name, offset: seq.GetOffset(), inset: seq.GetInset(), rc: false}
	return &seedSeq
}

//AddSeeds re-uses any seeds in the index, adding additional seeds as required to
//reach at least minSeeds for the sequence.
func (g *SeedIndex) AddSeeds(seq sequence.Sequence, minSeeds int, ranks []float64) {
	//count the number existing in the index and generate the gapped-seed sequence using those
	k := int(g.seedSize)

	//for speed, we duplicate Kmers() for a single pass
	var mask int
	for i := 0; i < k; i++ {
		mask = (mask << 2) | 3
	}

	count := seq.CountKmers(k, mask, g.seeds)

	if count < minSeeds {
		//find the n best seeds that will top us up to the required amount
		topN := make([]uint, minSeeds-count, minSeeds-count)
		topNValues := make([]float64, minSeeds-count, minSeeds-count)
		for i := 0; i < len(topNValues); i++ {
			topNValues[i] = math.MaxFloat64
		}
		//add additional seeds, not overlapping existing seeds
		kmer := seq.KmerAt(0, k)
		nextIndex := k
		for nextIndex < seq.Len()-k {
			//find the best k-mer in each k-length block. Resetting if a seed is found.
			reset := false
			bestValue := math.MaxFloat64
			var bestSeed int
			for i := 0; nextIndex < seq.Len() && i < k; i++ {
				kmer = seq.NextKmer(kmer, mask, nextIndex)
				nextIndex++
				if g.seeds[kmer] {
					reset = true
					break
				}
				value := ranks[kmer]
				if value < bestValue {
					bestValue = value
					bestSeed = kmer
				}
			}
			if !reset {
				n := 0
				for ; n < len(topNValues) && topNValues[n] > bestValue; n++ {
					if n > 0 {
						topNValues[n-1] = topNValues[n]
						topN[n-1] = topN[n]
					}
				}
				if n > 0 {
					topNValues[n-1] = bestValue
					topN[n-1] = uint(bestSeed)
				}

			}
			//step past the seed
			nextIndex += k
			if nextIndex < seq.Len()-k {
				kmer = seq.KmerAt(nextIndex, k)
			}
			nextIndex += k
		}
		// store the top seeds (synchronise here)
		g.lock.Lock()
		for _, seed := range topN {
			if !g.seeds[seed] {
				g.seeds[seed] = true
				g.size++
			}
			seed = ReverseComplement(seed, g.seedSize)
			if !g.seeds[seed] {
				g.seeds[seed] = true
				g.size++
			}
		}
		g.lock.Unlock()
	}
}

//NewAllSeedSequence creates a new SeedSequence adding every seed to the index
//Unlike NewSeedSequence this will not add the reversecomplement seeds to the index
func (g *SeedIndex) NewAllSeedSequence(seq sequence.Sequence) *SeedSequence {
	k := int(g.seedSize)
	var mask int
	for i := 0; i < k; i++ {
		mask = (mask << 2) | 3
	}
	segments := make([]int, 0, seq.Len()*2)
	prev := 0
	seed := seq.KmerAt(0, k) >> 2
	kmerIndex := 0
	g.lock.Lock()
	for i := k - 1; i < seq.Len(); i++ {
		seed = seq.NextKmer(seed, mask, i)
		if !g.seeds[seed] {
			g.seeds[seed] = true
			g.size++
		}
		segments = append(segments, kmerIndex-prev)
		segments = append(segments, seed)
		prev = kmerIndex + k
		kmerIndex++
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
		output <- index.NewSeedSequence(seq)
	}
	done <- true
}

//AddSeedsWorker processes all sequences, adding seeds to the provided index
func AddSeedsWorker(input <-chan sequence.Sequence, index *SeedIndex, numSeeds int, kmerValues []float64, done chan<- bool) {
	for seq := range input {
		index.AddSeeds(seq, numSeeds, kmerValues)
	}
	done <- true
}
