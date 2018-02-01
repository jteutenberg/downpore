package sequtil

import (
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/sequence"
	"sort"
)

func KmerOccurrences(seqs <-chan sequence.Sequence, k, numWorkers int) []uint64 {
	results := make(chan *[]uint64, numWorkers)
	for i := 0; i < numWorkers; i++ {
		go countWorker(seqs, k, results)
	}
	var counts []uint64
	for i := 0; i < numWorkers; i++ {
		if counts == nil {
			counts = *(<-results)
		} else {
			next := *(<-results)
			for j, c := range next {
				counts[j] += c
			}
		}
	}
	return counts
}

func countWorker(seqs <-chan sequence.Sequence, k int, done chan<- *[]uint64) {
	var mask int
	for i := 0; i < k; i++ {
		mask = (mask << 2) | 3
	}
	size := 1 << uint(k*2)
	counts := make([]uint64, size)
	for seq := range seqs {
		kmer := seq.KmerAt(0, k)
		counts[kmer]++
		for i := k; i < seq.Len(); i++ {
			kmer = seq.NextKmer(kmer, mask, i)
			counts[kmer]++
		}
	}
	done <- &counts
}

type sortableValues struct {
	values []uint64
	ids    []int
}

func (s *sortableValues) Len() int {
	return len(s.values)
}
func (s *sortableValues) Less(i, j int) bool {
	return s.values[i] < s.values[j]
}
func (s *sortableValues) Swap(i, j int) {
	s.ids[i], s.ids[j] = s.ids[j], s.ids[i]
	s.values[i], s.values[j] = s.values[j], s.values[i]
}

func TopOccurrences(counts []uint64, k uint, topN, bottomN int) (tops []int, bottoms []int) {
	//merge counts of forward and rc kmers
	ids := make([]int, len(counts))
	for i, c := range counts {
		ids[i] = i
		rc := seeds.ReverseComplement(uint(i), k)
		c += counts[rc]
		counts[i] = c
		counts[rc] = c

	}
	//sort
	vs := sortableValues{counts, ids}
	sort.Sort(&vs)

	//determine bottom from non-zero entries
	start := 0
	for start < len(counts) && vs.values[start] != 0 {
		start++
	}
	if start > len(counts)-bottomN {
		start = len(counts) - bottomN
	}
	//return slices of indices
	return vs.ids[start : start+bottomN], vs.ids[len(vs.ids)-topN:]
}
