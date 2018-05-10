package alignment

import (
	"fmt"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/util"
	"math"
)

//The initial offset is the position of each sequence's first kmer within the 32-kmer window
//Until the states have stepped this forward to position 16, no further window shifts are allowed
//There is also an embargo on landmarks until this process has completed.
//Basically, most pruning is turned off at the beginning while the sequences with various offsets sort themselves out.
var initialOffset = int32(8)

var debug = false

type Result struct {
	//the indices of each sequence that align to the last consensus base
	EndPositions []int32
}

//Aligner is an interface providing an Align method that is not thread-safe, it may (and should)
//re-use internal data-structures to minimise CPU/memory use
type Aligner interface {
	//GlobalConsensus finds the sequence that minimises the sum of alignment costs to it
	GlobalConsensus() (<-chan uint16, <-chan *QualityMetrics, <-chan *Result)
	//GlobalAlignment finds and returns the mutual alignment to the sequence that minimises alignment costs
	GlobalAlignment() (<-chan uint16, <-chan *QualityMetrics, <-chan []int)
	//GlobalAlignmentTo finds and returns the best mutual alignment to the given reference
	GlobalAlignmentTo([]uint16) (<-chan uint16, <-chan *QualityMetrics, <-chan []int)
	//ConsensusCost finds the cost of the best mutual alignment to the given reference
	ConsensusCost([]uint16) uint
}

type Measure interface {
	SetSequences(seqs [][]uint16, rc []bool)
	GetSequences() ([][]uint16, []bool)
	GetSequenceLen(int) int
	Distances(a uint16, sequence, start int, ds []uint16)
}

type QualityMetrics struct {
	ExactFraction  float64
	CostDelta      uint
	StateSpaceSize int
}

type simpleMeasure struct {
	seqs [][]uint16
	rcs  []bool
}

type dtw struct {
	maxWarp        int //the farthest divergence from the linear alignment (in bases)
	maxCost        uint16
	initialGapCost uint16
	costThreshold  uint16 //over a single sequence, for culling states
	measure        Measure
	full           bool //whether to enforce calling right to the end

	k                 int //k-mer length
	prevKmers         *util.IntSet
	kMask             uint16
	ds                []uint16
	landmarks         []*landmark
	expectedPositions []int32
	depth             int32
}

type landmark struct {
	k         uint16
	cost      uint
	seqs      []bool  //which sequences are restricted by this landmark
	positions []int32 //position in kmer sequence with the landmark kmer
}

type state struct {
	k            uint16     //the consensus k-mer
	positions    []int32    //which index we are up to within each sequence
	offsets      [][]uint16 //probability (cost) of each position (relative to lowest cost value)
	prev         *state     //prior state, one k-mer prior
	minCost      uint       //sum of lowest costs over all offsets
	votes        int        //number of exact matches for k
	spaceSize    int        //size of the state space when this was entered
	finished     bool       //if all seqs have reached their end (i.e. best pos is at len(seqs[i])-1)
	nextLandmark int        //which (current) landmark this is working towards
}

func NewFivemerMeasure() Measure {
	sm := simpleMeasure{}
	return &sm
}

func (m *simpleMeasure) SetSequences(seqs [][]uint16, rc []bool) {
	m.seqs = seqs
	m.rcs = rc
}
func (m *simpleMeasure) GetSequences() ([][]uint16, []bool) {
	return m.seqs, m.rcs
}
func (m *simpleMeasure) GetSequenceLen(index int) int {
	return len(m.seqs[index])
}
func (m *simpleMeasure) Distances(a uint16, sequence, start int, ds []uint16) {
	kmers := m.seqs[sequence]
	end := start + len(ds)
	if end > len(kmers) {
		for i := len(ds) + len(kmers) - end; i < len(ds); i++ {
			ds[i] = 1000
		}
		ds = ds[:len(ds)-end+len(kmers)]
		end = len(kmers)
	}
	for i := 0; i < len(ds); i++ {
		diff := kmers[start] ^ a
		//weighted sum of differences
		cost := (((diff >> 4) | (diff >> 5)) & 0x1) << 3 //centre: cost 8 for mismatch
		cost += (((diff >> 6) | (diff >> 7)) & 0x1) << 1 //one out: cost 2 for mismatch
		cost += (((diff >> 2) | (diff >> 3)) & 0x1) << 1
		cost += (((diff >> 1) | diff) & 0x1)    //edges: cost 1 for mismatch
		cost += (((diff >> 8) | diff>>9) & 0x1) //edges: cost 1 for mismatch
		ds[i] = cost
		start++
	}
}

func (s *state) printOffsets() {
	for j, offsets := range s.offsets {
		for _, off := range offsets {
			if off > 30000 {
				fmt.Print("- ")
			} else {
				fmt.Print(off, " ")
			}
		}
		fmt.Println("\n", s.positions[j])
	}
}

func NewDTWAligner(maxWarp int, initialGapCost uint, measure Measure, full bool, costThreshold uint, k int) Aligner {
	for maxWarp%8 != 0 {
		maxWarp++
	}
	ks := k
	var kMask uint16
	for k > 0 {
		kMask = (kMask << 2) | 0x3
		k--
	}
	d := dtw{maxWarp * 2, uint16(math.MaxUint16 / 2), uint16(initialGapCost), uint16(costThreshold), measure, full, ks, util.NewIntSet(), kMask, make([]uint16, maxWarp*2, maxWarp*2), make([]*landmark, 0, 10), nil, 0}
	return &d
}

func (lm *landmark) matchesLandmark(positions []int32) bool {
	count := 0
	c2 := 0
	for i, use := range lm.seqs {
		if use {
			c2++
			if positions[i] == lm.positions[i] {
				count++
			}
		}
	}
	return count >= c2/2
}
func (lm *landmark) matches(positions []int32) bool {
	count := 0
	c2 := 0
	for i, use := range lm.seqs {
		if use {
			c2++
			if positions[i] >= lm.positions[i]-6 && positions[i] <= lm.positions[i]+6 {
				count++
			}
		}
	}
	return count >= c2/2
}

//isPriorTo loosly checks for an ordering (nearby positions will be considered to be after the landmark)
func (lm *landmark) isPriorTo(positions []int32) bool {
	//landmark is prior if ALL relevant sequences are well before the given positions
	for i, use := range lm.seqs {
		if use && positions[i]-4 < lm.positions[i] {
			//this position might be before the landmark
			return false
		}
	}
	return true
}

//isPriorLandmarkTo is a strict ordering, assuming positions exactly match the landmark's kmer
func (lm *landmark) isPriorLandmarkTo(otherSeqs []bool, otherPositions []int32) bool {
	for i, use := range lm.seqs {
		if use && otherSeqs[i] && otherPositions[i] < lm.positions[i] {
			return false
		}
	}
	return true
}

func (lm *landmark) lockState(s *state, maxCost uint16) {
	centre := int32(len(s.offsets[0]) / 2)
	for j, p := range lm.positions {
		if !lm.seqs[j] || p < s.positions[j]-centre {
			continue
		}
		p -= s.positions[j] - centre
		offs := s.offsets[j]
		if int(p) < len(offs) {
			for n := 0; n < len(offs); n++ {
				if n != int(p) {
					offs[n] = maxCost
				}
			}
			s.minCost += uint(offs[p])
			offs[p] = 0
		}
	}
}

func (lm *landmark) cropState(s *state, maxCost uint16) {
	centre := int32(len(s.offsets[0]) / 2)
	for j, p := range lm.positions {
		if !lm.seqs[j] {
			continue
		}
		p -= s.positions[j] - centre
		offs := s.offsets[j]
		if p >= int32(len(offs)) || p < 0 {
			continue //this one didn't reach the landmark at all
		}
		for n := 0; n < int(p); n++ {
			offs[n] = maxCost
		}
		newMin := maxCost
		for n := int(p); n < len(offs); n++ {
			if offs[n] < newMin {
				newMin = offs[n]
			}
		}
		s.minCost += uint(newMin)
		for n := int(p); n < len(offs); n++ {
			offs[n] -= newMin
		}
	}
}

func fixDrift(s *state, bestPos, index int, maxCost uint16) {
	offs := s.offsets[index]
	centre := len(offs) / 2
	//have we drifted too far?
	drift := centre - bestPos
	if drift < -4 {
		//the best index is to the RHS: shuffle down
		for i := 0; i < len(offs)+drift; i++ {
			offs[i] = offs[i-drift]
		}
		s.positions[index] -= int32(drift)
		for i := len(offs) + drift; i < len(offs); i++ {
			offs[i] = maxCost
		}
	} else if drift > 4 {
		//shuffle up
		for i := len(offs) - 1; i >= drift; i-- {
			offs[i] = offs[i-drift]
		}
		s.positions[index] -= int32(drift)
		for i := 0; i < drift; i++ {
			offs[i] = maxCost
		}
	}
	//so each offset now contains the cheapest seq[offset]->kmer match available from all possible aligments that get there
}

func updateOffsetsAsm(ds, poffs, offsets []uint16, threshold uint16) uint16

func (d *dtw) prepareDistances(seq int, kmer uint16, pos, start, end int) (int, int) {
	centre := len(d.ds) / 2
	if start < 0 {
		start = 0
	}
	if end > len(d.ds) {
		end = len(d.ds)
	}
	seqStart := pos - centre + start
	if seqStart < 0 {
		start -= seqStart
		seqStart = 0
	}
	if start < 0 || start >= len(d.ds) || end < 0 || end > len(d.ds) || end <= start {
		fmt.Println("Over distance! ", start, end, " at pos ", pos)
	}
	d.measure.Distances(kmer, seq, seqStart, d.ds[start:end])
	for i := 0; i < start; i++ {
		d.ds[i] = 1256
	}
	for i := end; i < len(d.ds); i++ {
		d.ds[i] = 1256
	}
	//add distance from expected position
	exp := int(d.depth + d.expectedPositions[seq])
	for i := start; i < end; i++ {
		delta := (pos - centre + i) - exp
		if delta < -16 {
			d.ds[i] += uint16(-16 - delta)
		} else if delta > 16 {
			d.ds[i] += uint16(delta - 16)
		}
	}

	return start, end
}

func getBounds(values []uint16, maxValue uint16) (int, int) {
	start := 0
	for start < len(values) && values[start] >= maxValue {
		start++
	}
	end := len(values) - 1
	for end > 0 && values[end] >= maxValue {
		end--
	}
	end++
	start -= 2 //handle skips
	end++      //handle a stay
	return start, end
}

func getZeroPos(values []uint16, start, end int) int {
	for i := start; i < end; i++ {
		if values[i] == 0 {
			return i
		}
	}
	return len(values) / 2
}

func (d *dtw) updateCosts(s *state, prev *state, index int) (int, bool, uint16, bool) {
	//distances
	centre := len(s.offsets[index]) / 2
	pos := int(s.positions[index])
	//prev offsets are the step values. Compare against stay (prev offsets i+1)
	offsets := s.offsets[index]
	poffs := prev.offsets[index]
	//prepare the kmer distances
	start, end := getBounds(poffs, d.maxCost)
	if start < end && end >= 0 {
		//otherwise everything is over max cost -- don't need distances anyway
		start, end = d.prepareDistances(index, s.k, pos, start, end)
	}
	minCost := updateOffsetsAsm(d.ds, poffs, offsets, d.costThreshold)
	minPos := getZeroPos(offsets, start, end)
	exact := false
	for i := start; !exact && i < end; i++ {
		exact = d.ds[i] == 0 && offsets[i] < d.maxCost
	}
	if d.depth > initialOffset {
		fixDrift(s, minPos, index, d.maxCost)
	}
	return minPos, exact, uint16(minCost), pos+minPos-centre >= d.measure.GetSequenceLen(index)-1
}

func (s *state) traceBack(kmers chan uint16, costs chan *QualityMetrics) *state {
	cost := QualityMetrics{CostDelta: s.minCost}
	final := s
	if s.prev != nil {
		cost.CostDelta -= s.prev.minCost
		final = s.prev.traceBack(kmers, costs)
		s.prev = nil //help the GC
	}
	cost.ExactFraction = float64(s.votes) / float64(len(s.offsets))
	cost.StateSpaceSize = s.spaceSize
	kmers <- s.k
	costs <- &cost
	return final
}

func (s *state) traceBackFull(kmers chan uint16, costs chan *QualityMetrics, positions chan []int) *state {
	currentPos := make([]int, len(s.positions), len(s.positions))
	for i := 0; i < len(currentPos); i++ {
		bestPos := len(s.offsets[0]) - 1
		bestCost := s.offsets[i][bestPos]
		for j := bestPos - 1; j >= 0; j-- {
			if s.offsets[i][j] < bestCost {
				bestCost = s.offsets[i][j]
				bestPos = j
			}
		}
		//currentPos is position in the (whole) sequence
		currentPos[i] = int(s.positions[i]) + bestPos - len(s.offsets[i])/2
	}
	return s.traceBackFullAt(currentPos, kmers, costs, positions)
}

func (s *state) traceBackFullAt(currentPos []int, kmers chan uint16, costs chan *QualityMetrics, positions chan []int) *state {
	var cost QualityMetrics
	pos := make([]int, len(s.offsets), len(s.offsets))
	final := s
	for i, offsets := range s.offsets {
		bestCost := uint16(math.MaxUint16)
		bestPos := -1
		//find the index relating to the currentPos: we can't go beyond this
		latest := currentPos[i] - int(s.positions[i]) + len(offsets)/2 //current pos - position of index 0
		//find the cheapest base prior to the current position (actually, should be within 2 steps..)
		//fmt.Println("current pos:",currentPos[i],"state pos:",s.positions[i],"centre:",len(offsets)/2,"latest index:",latest)
		for j := latest; j >= latest-3 && j >= 0; j-- {
			if j >= len(offsets) {
				continue
			}
			c := offsets[j]
			if c < bestCost {
				bestCost = c
				bestPos = j
			}
		}
		bestPos += int(s.positions[i]) - len(offsets)/2 //turn the index into a position
		pos[i] = bestPos
	}
	cost.CostDelta = s.minCost
	if s.prev != nil {
		cost.CostDelta -= s.prev.minCost
		final = s.prev.traceBackFullAt(pos, kmers, costs, positions)
		s.prev = nil //help the GC
	}
	cost.ExactFraction = float64(s.votes) / float64(len(s.offsets))
	cost.StateSpaceSize = s.spaceSize
	kmers <- s.k
	costs <- &cost
	positions <- pos
	return final
}

func (d *dtw) nextState(current []*state, next *[]*state, nextK uint16) bool {
	d.depth++
	s := current[0]
	if s.finished {
		*next = append(*next, s)
		return true
	}
	successor := state{nextK, make([]int32, len(s.positions), len(s.positions)), make([][]uint16, len(s.offsets), len(s.offsets)), s, s.minCost, 1, 1, false, s.nextLandmark}
	for j := 0; j < len(successor.offsets); j++ {
		successor.offsets[j] = make([]uint16, len(s.offsets[j]), len(s.offsets[j]))
	}
	if d.full {
		successor.finished = true
	}
	var tailGap uint
	for j, p := range s.positions {
		successor.positions[j] = p + 1
		_, _, cost, finished := d.updateCosts(&successor, s, j)
		successor.minCost += uint(cost)
		if !finished {
			tailGap += uint(d.measure.GetSequenceLen(j) - 1 - int(successor.positions[j]))
		}
		/*if debug {
			for _, off := range successor.offsets[j] {
				if off > 30000 {
					fmt.Print("- ")
				} else {
					fmt.Print(off," ")
				}
			}
			fmt.Println("\n",successor.positions[j],ex)
		}*/
		if d.full {
			successor.finished = successor.finished && finished
		} else {
			successor.finished = successor.finished || finished
		}
	}
	if successor.finished {
		successor.minCost += tailGap * uint(d.initialGapCost)
	}
	//if debug{
	//	fmt.Println(sequence.KmerString(int(successor.k),int(d.k)))
	//}
	*next = append(*next, &successor)
	return successor.finished
}

func (d *dtw) nextStates(current []*state, next *[]*state) bool {
	d.depth++
	d.prevKmers.Clear()
	minFinishedCost := uint(math.MaxUint32)
	allFinished := true
	highCounts := 0
	minLandmark := 100000
	landmarkAdded := false
	lowestCost := uint(math.MaxUint32)
	for _, s := range current {
		if (len(d.landmarks) == 0 || s.nextLandmark == len(d.landmarks)) && s.minCost < lowestCost {
			lowestCost = s.minCost
		}
		if s.finished && s.minCost < minFinishedCost {
			minFinishedCost = s.minCost
		}
	}
	maxDelta := uint(len(current[0].offsets)/2) * uint(d.costThreshold)
	lowestCost += maxDelta
	for m, s := range current {
		if s.finished {
			//keep hold of all the finished states. Let the remainder look a bit further ahead.
			if minFinishedCost >= s.minCost {
				*next = append(*next, s)
			}
			if debug {
				fmt.Println("Finished", sequence.KmerString(int(s.k), int(d.k)), "keep =", (minFinishedCost >= s.minCost))
			}
			continue
		}
		if s.minCost > lowestCost {
			if debug {
				fmt.Println("Discarding high cost:", s.minCost, "vs", lowestCost, ":", sequence.KmerString(int(s.k), int(d.k)))
			}
			continue
		}
		shifted := (s.k << 2) & d.kMask
		iShift := uint(shifted)
		update := d.prevKmers.Contains(iShift)
		added := false
		if debug {
			fmt.Println("SUCCESSORS OF", sequence.KmerString(int(s.k), int(d.k)), "at", s.nextLandmark, "/", len(d.landmarks))
		}
		//four possible steps of k-mer
		for i := uint16(0); i < 4; i++ {
			nextK := shifted | i
			if nextK == s.k { //homopolymer
				continue
			}
			successor := state{nextK, make([]int32, len(s.positions), len(s.positions)), make([][]uint16, len(s.offsets), len(s.offsets)), s, s.minCost, 0, 0, false, s.nextLandmark}
			for j := 0; j < len(successor.offsets); j++ {
				successor.offsets[j] = make([]uint16, len(s.offsets[j]), len(s.offsets[j]))
			}

			votes := 0 //number of exact matches
			lastVoted := -1
			lastVotedIndex := -1
			var extraCost uint
			if d.full {
				successor.finished = true
			}
			for j, p := range s.positions {
				successor.positions[j] = p + 1
				minIndex, exactMatch, cost, finished := d.updateCosts(&successor, s, j)
				if exactMatch {
					votes++
					lastVoted = j
					lastVotedIndex = minIndex
				}
				extraCost += uint(cost)
				if d.full {
					successor.finished = successor.finished && finished
				} else {
					successor.finished = successor.finished || finished
				}
			}
			//successor.minCost += (2*extraCost) / uint(votes+1)
			successor.minCost += extraCost
			successor.votes = votes
			if successor.finished {
				if minFinishedCost > successor.minCost {
					minFinishedCost = successor.minCost
				}
			}
			if debug {
				fmt.Println(i, "child present in", votes, "and min cost up by", extraCost, "from", s.minCost, "=", successor.minCost)
			}
			if votes == 0 { //not present in any sequence at any relevant position
				continue
			}
			if votes == 1 {
				//remove all from offset except the exact match
				successor.minCost += uint(successor.offsets[lastVoted][lastVotedIndex])
				for n, _ := range successor.offsets[lastVoted] {
					if n != lastVotedIndex {
						successor.offsets[lastVoted][n] = d.maxCost
					} else {
						successor.offsets[lastVoted][n] = 0
					}
				}
			}
			//test this sequence against existing landmarks. Not allowed to skip over any. If this is a landmark of its own, handle later.
			if successor.nextLandmark < len(d.landmarks) {
				lm := d.landmarks[successor.nextLandmark]
				if successor.minCost > lm.cost {
					//already more expensive than its upcoming alternative at the landmark
					if debug {
						fmt.Println("Cost past next landmark", successor.nextLandmark, "/", len(d.landmarks), ":", successor.minCost, lm.cost)
					}
					continue
				}
				if nextK == lm.k && lm.matches(successor.positions) {
					if debug {
						fmt.Println("Reached my next landmark. Crop=", votes <= len(successor.offsets)/2)
					}
					//we're at the landmark, but only continue if we are a better alternative??
					if votes <= len(successor.offsets)/2 {
						lm.cropState(&successor, d.maxCost) //make sure we can't jump back before it again
					}
					successor.nextLandmark++
				} else if lm.isPriorTo(successor.positions) {
					continue //discard this landmark violator!
				}
			}
			if !successor.finished && d.depth > initialOffset && votes > len(s.offsets)/2 {
				lmPositions := make([]int32, len(s.offsets), len(s.offsets))
				lmSeq := make([]bool, len(s.offsets), len(s.offsets))
				lmCost := successor.minCost
				seqs, _ := d.measure.GetSequences()
				for j, pos := range successor.positions {
					seq := seqs[j]
					seqLen := int32(len(seq))
					offs := successor.offsets[j]
					if pos < seqLen && seq[pos] == nextK {
						lmPositions[j] = pos
						off := offs[len(offs)/2]
						lmSeq[j] = off < d.maxCost
						if off < d.maxCost {
							lmCost += uint(off)
						}
					} else {
						for k := int32(1); k < 16; k++ {
							if pos+k < seqLen && seq[pos+k] == nextK {
								lmPositions[j] = pos + k
								off := offs[len(offs)/2+int(k)]
								lmSeq[j] = off < d.maxCost
								if off < d.maxCost {
									lmCost += uint(off)
								}
								break
							}
							if pos-k >= 0 && pos-k < seqLen && seq[pos-k] == nextK {
								lmPositions[j] = pos - k
								off := offs[len(offs)/2-int(k)]
								lmSeq[j] = off < d.maxCost
								if off < d.maxCost {
									lmCost += uint(off)
								}
								break
							}
						}
					}
				}
				if debug {
					fmt.Println("LANDMARK voted", votes, sequence.KmerString(int(nextK), int(d.k)), "at", lmPositions, lmCost)
				}
				//this landmark needs to be before my next one
				if successor.nextLandmark < len(d.landmarks) && d.landmarks[successor.nextLandmark].isPriorLandmarkTo(lmSeq, lmPositions) {
					continue
				}
				var mark *landmark
				deleteAfter := true
				updatedLandmark := false //have I reduced the cost of an existing landmark?
				skippedLandmark := false //have I gone past a required landmark?
				//test for existing landmark, include one prior
				if len(d.landmarks) > 0 {
					j := 0
					if successor.nextLandmark > 0 {
						j = successor.nextLandmark - 1
					}
					for ; j < len(d.landmarks); j++ {
						lm := d.landmarks[j]
						if lm.k == nextK && lm.matchesLandmark(lmPositions) {
							skippedLandmark = skippedLandmark || successor.nextLandmark < j
							mark = lm
							if debug {
								fmt.Println("matches an existing landmark at", j, "my next was", successor.nextLandmark, ". cost:", lmCost, "vs", lm.cost, sequence.KmerString(int(nextK), int(d.k)), "and", sequence.KmerString(int(lm.k), int(d.k)))
							}
							//delete me or update the landmark? Only if I've passed all other landmarks
							if !skippedLandmark && lm.cost > lmCost {
								lm.cost = lmCost
								lm.positions = lmPositions
								lm.seqs = lmSeq
								lm.lockState(&successor, d.maxCost)
								//all later landmarks are invalid now
								d.landmarks = d.landmarks[:j+1]
								updatedLandmark = true
							} else {
								successor.nextLandmark = j + 1
								deleteAfter = false
							}
							break
						}
					}
				}
				if skippedLandmark {
					continue
				}
				if mark == nil { //a brand new landmark
					lm := landmark{k: nextK, cost: lmCost, seqs: lmSeq, positions: lmPositions}
					mark = &lm
					//any later landmarks are irrelevant now
					newLen := len(d.landmarks)
					for newLen > 0 && mark.isPriorLandmarkTo(d.landmarks[newLen-1].seqs, d.landmarks[newLen-1].positions) {
						newLen--
					}
					d.landmarks = d.landmarks[:newLen]
					d.landmarks = append(d.landmarks, mark)
					successor.nextLandmark = len(d.landmarks)
					if debug {
						fmt.Println("adding new landmark", lmPositions)
					}
					//collapse self to only the landmark position
					mark.lockState(&successor, d.maxCost)
					landmarkAdded = true
				}
				if deleteAfter {
					//all states after this have not been through this landmark. Purge!
					for j := len(*next) - 1; j >= 0; j-- {
						n := (*next)[j]
						//if this is my previously achieved landmark but with an updated cost.. I should be removed
						if (updatedLandmark && n.nextLandmark >= len(d.landmarks)-1) || mark.isPriorTo(n.positions) || n.minCost > mark.cost {
							(*next)[j] = (*next)[len(*next)-1]
							*next = (*next)[:len(*next)-1]
						} else if n.nextLandmark > len(d.landmarks)-1 {
							//otherwise, they have to achieve this landmark at some point
							n.nextLandmark = len(d.landmarks) - 1
						}
					}
					//any upcoming states must be prior to this too
					for j := len(current) - 1; j >= m+1; j-- {
						if current[j].nextLandmark >= len(d.landmarks)-1 && current[j].k == mark.k && mark.matches(current[j].positions) {
							current[j].nextLandmark = len(d.landmarks) //achieved the new landmark just now!
							mark.cropState(current[j], d.maxCost)      //make sure we can't jump back before it again
						} else if (updatedLandmark && current[j].nextLandmark >= len(d.landmarks)-1) || current[j].minCost > mark.cost || mark.isPriorTo(current[j].positions) {
							current[j] = current[len(current)-1]
							current = current[:len(current)-1]
						} else if current[j].nextLandmark >= len(d.landmarks) {
							current[j].nextLandmark = len(d.landmarks) - 1
						}
					}
				} else {
					continue //this matched another landmark and was more expensive
				}
			}
			if minFinishedCost >= successor.minCost {
				added = true
				if update {
					found := false
					keep := false
					//find the existing target next state and consider replacing it
					for j, other := range *next {
						if other.k == nextK {
							found = true
							//is the other one more expensive (and I'm not waiting on its landmark)
							if other.minCost >= successor.minCost && other.nextLandmark <= successor.nextLandmark {
								//do the substitution
								(*next)[j] = &successor
							} else {
								//keep me if I'm ahead of its landmark -- it is probably an invalid state
								keep = keep || other.nextLandmark < successor.nextLandmark
							}
						}
					}
					if !found || keep {
						//this will occur if the last version had > minFinishedCost
						allFinished = false
						*next = append(*next, &successor)
						if votes > len(successor.offsets)/4 {
							highCounts++
						}
						if successor.nextLandmark < minLandmark {
							minLandmark = successor.nextLandmark
						}
					}
				} else {
					allFinished = false
					*next = append(*next, &successor)
					if votes > len(successor.offsets)/4 {
						highCounts++
					}
					if successor.nextLandmark < minLandmark {
						minLandmark = successor.nextLandmark
					}
				}

			}
		}
		//if any one of the four has 2x the number of votes of another (and >50% possible votes), then remove the other states
		if !update && added {
			d.prevKmers.Add(iShift)
		}
	}
	// perform vote pruning?
	if debug {
		fmt.Println("have ", len(*next), "states with", highCounts, "high vote states. min landmark is", minLandmark, "/", len(d.landmarks))
	}
	if d.depth > initialOffset && len(*next) > 20 && highCounts > 5 {
		for j := len(*next) - 1; j >= 0; j-- {
			n := (*next)[j]
			if n.nextLandmark == minLandmark && (n.votes < 2 || n.votes < len(n.offsets)/8) {
				(*next)[j] = (*next)[len(*next)-1]
				*next = (*next)[:len(*next)-1]
			}
		}
	}
	if landmarkAdded {
		d.updateExpectedPositions()
	}
	sSize := len(*next)
	for _, s := range *next {
		s.spaceSize = sSize
	}
	return allFinished
}

func (d *dtw) newState(k uint16) *state {
	seqs, _ := d.measure.GetSequences()
	s := state{k, make([]int32, len(seqs), len(seqs)), make([][]uint16, len(seqs), len(seqs)), nil, 0, 0, 0, false, 0}
	//setup the offsets
	for i, seq := range seqs {
		s.positions[i] = int32(initialOffset)
		s.offsets[i] = make([]uint16, d.maxWarp, d.maxWarp)
		if seq[0] != k {
			s.offsets[i][initialOffset] = d.initialGapCost
		} else {
			s.offsets[i][initialOffset] = 0
		}
		for j := int(initialOffset + 1); j < len(s.offsets[i]); j++ {
			s.offsets[i][j] = d.initialGapCost
		}
		for j := 0; j < int(initialOffset); j++ {
			s.offsets[i][j] = d.maxCost
		}
	}
	return &s
}

func (d *dtw) firstStates() []*state {
	states := make([]*state, 0, 100)
	d.prevKmers.Clear()
	seqs, _ := d.measure.GetSequences()
	for _, seq := range seqs {
		d.prevKmers.Add(uint(seq[0]))
	}
	for ok, nextK := d.prevKmers.GetFirstID(); ok; ok, nextK = d.prevKmers.GetNextID(nextK) {
		k := uint16(nextK)
		states = append(states, d.newState(k))
	}
	for _, s := range states {
		s.spaceSize = len(states)
	}
	return states
}

//writeBestPositions overwrites the state's position slice with pos+offset values
func (s *state) writeBestPositions() {
	for i, offs := range s.offsets {
		bp := -1
		bestCost := uint16(math.MaxUint16)
		for j, off := range offs {
			if off < bestCost {
				bestCost = off
				bp = j
			}
		}
		//re-use the positions slice
		s.positions[i] += int32(bp - len(offs)/2)
	}
}

//updateExpectedPositions determines where in the sequence we expect to be, assuming
//balanced inserts and deletes from the last landmark
func (d *dtw) updateExpectedPositions() {
	//assume we have a new landmark
	lm := d.landmarks[len(d.landmarks)-1]
	for i, use := range lm.seqs {
		if use {
			d.expectedPositions[i] = lm.positions[i] - d.depth
		}
	}
}

//GlobalConsensus returns a consensus sequence as k-mer, cost pairs
func (d *dtw) GlobalConsensus() (<-chan uint16, <-chan *QualityMetrics, <-chan *Result) {
	d.depth = 0
	//"states" will form a tree of consensus calls
	//each state contains a small band of a probability distribution over each sequence's position
	basesOut := make(chan uint16, 20)
	costsOut := make(chan *QualityMetrics, 20)
	finalResult := make(chan *Result, 1)
	seqs, _ := d.measure.GetSequences()
	d.expectedPositions = make([]int32, len(seqs), len(seqs))

	states := d.firstStates()

	next := make([]*state, 0, 100)
	finished := false

	go func() {
		for !finished {
			if debug {
				fmt.Println("COUNT ", d.depth)
			}
			finished = d.nextStates(states, &next)
			if debug {
				fmt.Println()
				fmt.Println("Length", d.depth, "have", len(next), "next states....")
			}
			//fmt.Println(d.expectedPositions)
			//for _, s := range next {
			//	fmt.Println(sequence.KmerString(int(s.k),int(d.k)),"(",s.nextLandmark,")\n",s.positions)
			//}
			/*for _, s := range next {
				fmt.Print(sequence.KmerString(int(s.k),int(d.k)),"(",s.nextLandmark,"),")
			}
			fmt.Println()*/
			//for _, lm := range d.landmarks {
			//	fmt.Print(sequence.KmerString(int(lm.k),int(d.k)),"->")
			//}
			//fmt.Println()

			if !finished && len(next) == 1 && next[0].prev != nil {
				next[0].prev.traceBack(basesOut, costsOut)
				next[0].prev = nil
			}
			if len(next) == 0 {
				break
			}
			next, states = states[:0], next
		}
		if len(states) > 0 {
			best := -1
			bestCost := uint(math.MaxUint32)
			for i, s := range states {
				if s.minCost < bestCost {
					best = i
					bestCost = s.minCost
				}
			}
			firstState := states[best].traceBack(basesOut, costsOut)
			//find best positions for each sequence
			states[best].writeBestPositions()
			firstState.writeBestPositions()
			//fmt.Println(states[best].positions)
			result := Result{EndPositions: states[best].positions}
			finalResult <- &result
		}

		close(basesOut)
		close(costsOut)
		close(finalResult)
	}()
	return basesOut, costsOut, finalResult
}

func (d *dtw) GlobalAlignment() (<-chan uint16, <-chan *QualityMetrics, <-chan []int) {
	d.depth = 0
	basesOut := make(chan uint16, 20)
	costsOut := make(chan *QualityMetrics, 20)
	posOut := make(chan []int, 20)
	seqs, _ := d.measure.GetSequences()
	d.expectedPositions = make([]int32, len(seqs), len(seqs))

	states := d.firstStates()
	next := make([]*state, 0, 100)
	finished := false

	go func() {
		for !finished {
			finished = d.nextStates(states, &next)
			if !finished && len(next) == 1 && next[0].prev != nil {
				next[0].prev.traceBackFull(basesOut, costsOut, posOut)
				next[0].prev = nil
			}
			if len(next) == 0 {
				break
			}
			next, states = states[:0], next
		}
		if len(states) > 0 {
			best := -1
			bestCost := uint(math.MaxUint32)
			for i, s := range states {
				if s.minCost < bestCost {
					best = i
					bestCost = s.minCost
				}
			}
			states[best].traceBackFull(basesOut, costsOut, posOut)
			//fmt.Println(states[best].positions)
		}
		close(posOut)
		close(basesOut)
		close(costsOut)
	}()
	return basesOut, costsOut, posOut
}

func (d *dtw) GlobalAlignmentTo(reference []uint16) (<-chan uint16, <-chan *QualityMetrics, <-chan []int) {
	d.depth = 0
	basesOut := make(chan uint16, 20)
	costsOut := make(chan *QualityMetrics, 20)
	posOut := make(chan []int, 20)
	seqs, _ := d.measure.GetSequences()
	d.expectedPositions = make([]int32, len(seqs), len(seqs))

	states := make([]*state, 1, 2)
	states[0] = d.newState(reference[0])
	states[0].spaceSize = 1
	next := make([]*state, 0, 2)
	finished := false

	go func() {
		for i := 1; i < len(reference) && !finished; i++ {
			finished = d.nextState(states, &next, reference[i])
			next, states = states[:0], next
		}
		states[0].traceBackFull(basesOut, costsOut, posOut)
		close(basesOut)
		close(posOut)
		close(costsOut)
	}()
	return basesOut, costsOut, posOut
}

func (d *dtw) ConsensusCost(reference []uint16) uint {
	d.depth = 0
	states := make([]*state, 1, 2)
	states[0] = d.newState(reference[0])
	states[0].spaceSize = 1
	next := make([]*state, 0, 2)
	seqs, _ := d.measure.GetSequences()
	d.expectedPositions = make([]int32, len(seqs), len(seqs))
	finished := false
	for i := 1; i < len(reference) && !finished; i++ {
		finished = d.nextState(states, &next, reference[i])
		next, states = states[:0], next
	}
	return uint(states[0].minCost)
}
