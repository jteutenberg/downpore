package seeds

import (
	"fmt"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/util"
	"sync"
)

//SeedIndex is a set of reads in gapped-seed format ready for querying
type SeedIndex struct {
	seedSize     uint
	kmers        []bool          //which k-mers have seeds in this index
	sequences    []*SeedSequence //list of all sequences
	sequenceSets []*util.IntSet  //list of seed -> set of sequences (indices)
	seedSets     []*util.IntSet  //list of sequence -> set of seeds
	kmerMap      []int32         //maps kmer->seed
	seedMap      []int           //maps seed->kmer
	size         int
	lock         *sync.Mutex
}

func NewSeedIndex(k uint) *SeedIndex {
	size := 1
	for j := k; j > 0; j-- {
		size *= 4
	}
	var lock sync.Mutex
	sg := SeedIndex{kmers: make([]bool, size), kmerMap: make([]int32, size), seedSize: k, sequences: make([]*SeedSequence, 0, 100000), sequenceSets: make([]*util.IntSet, 0, 100000), seedSets: make([]*util.IntSet, 0, 100000), size: 0, lock: &lock}
	return &sg
}

func (g *SeedIndex) NewSeedSequence(seq sequence.Sequence) *SeedSequence {
	k := int(g.seedSize)
	//for speed, we duplicate Kmers() for a single pass
	var mask int
	for i := 0; i < k; i++ {
		mask = (mask << 2) | 3
	}
	count := seq.CountKmers(seq.Len(), k, mask, g.kmers)
	segments := make([]int, count*2+1)
	seq.WriteSegments(segments, k, mask, g.kmers)
	//translate the k-mers into seeds
	for i := 1; i < len(segments); i += 2 {
		segments[i] = int(g.kmerMap[segments[i]])
	}
	name := seq.GetName()
	seedSeq := SeedSequence{segments: segments, length: seq.Len(), id: seq.GetID(), name: &name, offset: seq.GetOffset(), inset: seq.GetInset(), rc: false}
	return &seedSeq
}

func (g *SeedIndex) SeedString(seed int) string {
	return sequence.KmerString(g.seedMap[seed], int(g.seedSize))
}

func (g *SeedIndex) SeedCount(seed int) uint {
	return g.sequenceSets[seed].Size()
}

//AddSeeds re-uses any seeds in the index, adding additional seeds as required to
//reach at least minSeeds for the sequence.
func (g *SeedIndex) AddSeeds(seq sequence.Sequence, minSeeds int, kmerRanks []float64) {
	//count the number existing in the index and generate the gapped-seed sequence using those
	k := int(g.seedSize)

	//for speed, we duplicate Kmers() for a single pass
	var mask int
	for i := 0; i < k; i++ {
		mask = (mask << 2) | 3
	}

	count := seq.CountKmers(minSeeds, k, mask, g.kmers)
	q := seq.Quality()
	count = 0 //a full set of seeds for this sequence

	if count < minSeeds {
		//find the n best seeds that will top us up to the required amount
		topN := make([]uint, minSeeds-count, minSeeds-count)
		topNValues := make([]float64, minSeeds-count, minSeeds-count)
		for i := 0; i < len(topNValues); i++ {
			topNValues[i] = 0
		}
		//add additional seeds, not overlapping existing seeds
		kmer := seq.KmerAt(0, k)
		nextIndex := k
		for nextIndex < seq.Len()-k {
			//find the best k-mer in each k-length block. Resetting if a seed is found.
			reset := false
			bestValue := 0.0
			var bestSeed int
			for i := 0; nextIndex < seq.Len() && i < k; i++ {
				kmer = seq.NextKmer(kmer, mask, nextIndex)
				nextIndex++
				if g.kmers[kmer] {
					reset = true
					break
				}
				value := kmerRanks[kmer]
				if q != nil {
					value *= float64(q[nextIndex-k/2])
				}
				if value > bestValue {
					bestValue = value
					bestSeed = kmer
				}
			}
			if !reset {
				n := 0
				//position 0 is the bottom spot. Shuffle down to make room if needed.
				for ; n < len(topNValues) && topNValues[n] < bestValue; n++ {
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
		for _, kmer := range topN {
			if !g.kmers[kmer] {
				g.kmers[kmer] = true
				g.kmerMap[kmer] = int32(g.size)
				for len(g.sequenceSets) <= g.size {
					g.sequenceSets = append(g.sequenceSets, util.NewIntSet())
					g.seedMap = append(g.seedMap, -1)
				}
				g.seedMap[g.size] = int(kmer)
				g.size++
			}
			kmer = ReverseComplement(kmer, g.seedSize)
			if !g.kmers[kmer] {
				g.kmers[kmer] = true
				g.kmerMap[kmer] = int32(g.size)
				for len(g.sequenceSets) <= g.size {
					g.sequenceSets = append(g.sequenceSets, util.NewIntSet())
					g.seedMap = append(g.seedMap, -1)
				}
				g.seedMap[g.size] = int(kmer)
				g.size++
			}
		}
		g.lock.Unlock()
	}
}

//AddSingleSeeds adds one seet for every seedRate bases if no seeds already exists within them.
//The seed with minimal rank is always selected.
func (g *SeedIndex) AddSingleSeeds(seq sequence.Sequence, seedRate int, ranks []float64) {
	//count the number existing in the index and generate the gapped-seed sequence using those
	k := int(g.seedSize)

	var mask int
	for i := 0; i < k; i++ {
		mask = (mask << 2) | 3
	}

	//TODO: feed each subsequence to a goroutine
	for i := 0; i < seq.Len()-seedRate; i += seedRate {
		count := seq.CountKmersBetween(i, i+seedRate, 1, k, mask, g.kmers)
		if count == 0 {
			end := i + seedRate
			//find the  best seed
			kmer := seq.KmerAt(i, k)
			bestValue := ranks[kmer]
			bestKmer := kmer
			for j := i + k; j < end; j++ {
				kmer = seq.NextKmer(kmer, mask, j)
				value := ranks[kmer]
				if value > bestValue {
					bestValue = value
					bestKmer = kmer
				}
			}
			g.lock.Lock()
			if !g.kmers[bestKmer] { //in case of race condition
				g.kmers[bestKmer] = true
				g.kmerMap[bestKmer] = int32(g.size)
				for len(g.sequenceSets) <= g.size {
					g.sequenceSets = append(g.sequenceSets, util.NewIntSet())
					g.seedMap = append(g.seedMap, -1)
				}
				g.seedMap[g.size] = bestKmer
				g.size++
			}
			g.lock.Unlock()
		}
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
	kmer := seq.KmerAt(0, k) >> 2
	kmerIndex := 0
	g.lock.Lock()
	for i := k - 1; i < seq.Len(); i++ {
		kmer = seq.NextKmer(kmer, mask, i)
		if !g.kmers[kmer] {
			g.kmers[kmer] = true
			g.kmerMap[kmer] = int32(g.size)
			for len(g.sequenceSets) <= g.size {
				g.sequenceSets = append(g.sequenceSets, util.NewIntSet())
				g.seedMap = append(g.seedMap, -1)
			}
			g.seedMap[g.size] = kmer
			g.size++
		}
		segments = append(segments, kmerIndex-prev)
		segments = append(segments, int(g.kmerMap[kmer]))
		prev = kmerIndex + k
		kmerIndex++
	}
	g.lock.Unlock()
	segments = append(segments, 0)
	name := seq.GetName()
	seedSeq := SeedSequence{segments: segments, length: seq.Len(), id: seq.GetID(), name: &name, offset: seq.GetOffset(), inset: seq.GetInset(), rc: false}
	return &seedSeq
}

func (g *SeedIndex) Size() int {
	return g.size
}

func (g *SeedIndex) GetSeedSequence(index uint) *SeedSequence {
	return g.sequences[index]
}

func (g *SeedIndex) GetSeedsFromKmers(kmers []uint16, seedSet *util.IntSet) {
	for _, k := range kmers {
		if g.kmers[k] {
			seedSet.Add(uint(g.kmerMap[k]))
		}
	}
}

//GetSeedSet returns a set of seeds for the given sequence
func (g *SeedIndex) GetSeedSet(index uint) *util.IntSet {
	return g.seedSets[index]
}

func (g *SeedIndex) GetSeedLength() uint {
	return g.seedSize
}

func (g *SeedIndex) GetNumSequences() uint {
	return uint(len(g.sequences))
}

/*func (g *SeedIndex) Contains(kmer uint) bool {
	return g.kmers[kmer]
}*/

func (g *SeedIndex) AddSequence(seq *SeedSequence) {
	maxSeed := 0
	for i := 1; i < len(seq.segments); i += 2 {
		seed := seq.segments[i]
		//g.sequenceSets[seed].Add(index)
		if seed > maxSeed {
			maxSeed = seed
		}
	}
	seedSet := util.NewIntSetCapacity(maxSeed + 1)
	for i := 1; i < len(seq.segments); i += 2 {
		seed := seq.segments[i]
		seedSet.Add(uint(seed))
	}
	g.lock.Lock()
	g.sequences = append(g.sequences, seq)
	g.seedSets = append(g.seedSets, seedSet)
	g.lock.Unlock()
}

func (g *SeedIndex) IndexSequences(numWorkers int) {
	done := make(chan bool, numWorkers)
	start := 0
	split := len(g.sequenceSets) / numWorkers
	for i := 0; i < numWorkers-1; i++ {
		end := start + split - 1
		go g.index(start, end, done)
		start = end + 1
	}
	go g.index(start, len(g.sequenceSets), done)
	for i := 0; i < numWorkers; i++ {
		<-done
	}
}

func (g *SeedIndex) RemoveSequences() {
	g.sequences = g.sequences[:0]
	g.seedSets = g.seedSets[:0]
	for _, seqs := range g.sequenceSets {
		if seqs != nil {
			seqs.Clear()
		}
	}
}

func (g *SeedIndex) Destroy() {
	g.sequences = nil
	g.sequenceSets = nil
	g.seedSets = nil
	g.kmers = nil
	g.seedMap = nil
	g.kmerMap = nil
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
	return util.GetSharedIDs(allSeedSets, minCount, true)
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

func (g *SeedIndex) index(minSeed, maxSeed int, done chan<- bool) {
	for i := len(g.sequences) - 1; i >= 0; i-- {
		ind := uint(i)
		s := g.sequences[ind]
		for j := 1; j < len(s.segments); j += 2 {
			seed := s.segments[j]
			if seed >= minSeed && seed <= maxSeed {
				g.sequenceSets[seed].Add(ind)
			}
		}
	}
	done <- true
}
