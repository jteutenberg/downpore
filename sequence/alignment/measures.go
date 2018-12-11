package alignment

import (
	//"fmt"
	// "github.com/jteutenberg/downpore/sequence"
 )

type simpleMeasure struct {
	seqs [][]uint16
	rcs  []bool
	k int
}

func NewFivemerMeasure() Measure {
	sm := simpleMeasure{k:5}
	return &sm
}
func NewThreemerMeasure() Measure {
	sm := simpleMeasure{k:3}
	return &sm
}
func NewBaseMeasure() Measure {
	sm := simpleMeasure{k:1}
	return &sm
}
func NewFourmerMeasure() Measure {
	sm := simpleMeasure{k:4}
	return &sm
}
func NewSixmerMeasure() Measure {
	sm := simpleMeasure{k:6}
	return &sm
}

func (m *simpleMeasure) SetSequences(seqs [][]uint16, rc []bool) {
	m.seqs = seqs
	m.rcs = rc
}
func (m *simpleMeasure) GetSequences() ([][]uint16, []bool) {
	return m.seqs,m.rcs
}
func (m *simpleMeasure) GetSequenceLen(index int) int {
	return len(m.seqs[index])
}
func (m *simpleMeasure) Distances(a uint16, sequence, start int, ds []uint16) {
	kmers := m.seqs[sequence]
	end := start + len(ds)
	if end > len(kmers) {
		f := len(ds) + len(kmers) - end
		if f < 0 {
			f = 0
		}
		for i := f; i < len(ds); i++ {
			ds[i] = 14 //a bit high for smaller k (lower max cost)
		}
		ds = ds[:f]
	}
	// the dirty work
	if m.k == 5 {
		for i := 0; i < len(ds); i++ {
			diff := kmers[start] ^ a
			//weighted sum of differences
			cost := (((diff >> 4) | (diff >> 5)) & 0x1) << 3 //centre: cost 8 for mismatch
			cost += (((diff >> 6) | (diff >> 7)) & 0x1) << 1 //one out: cost 2 for mismatch
			cost += (((diff >> 2) | (diff >> 3)) & 0x1) << 1
			cost += (((diff >> 1) | diff) & 0x1)    //edges: cost 1 for mismatch
			cost += (((diff >> 8) | diff>>9) & 0x1)
			ds[i] = cost
			start++
		}
	} else if m.k == 4 {
		for i := 0; i < len(ds); i++ {
			diff := kmers[start] ^ a
			cost := (((diff >> 4) | (diff >> 5)) & 0x1) << 2 //centres: cost 4 for mismatch
			cost += (((diff >> 2) | (diff >> 3)) & 0x1) << 2 //centres: cost 4 for mismatch
			cost += (((diff >> 6) | (diff >> 7)) & 0x1) << 1 //one out: cost 2 for mismatch
			cost += (((diff >> 1) | diff) & 0x1) << 1
			ds[i] = cost
			start++
		}

	} else if m.k == 3 {
		for i := 0; i < len(ds); i++ {
			diff := kmers[start] ^ a
			cost := (((diff >> 2) | (diff >> 3)) & 0x1) << 3 //centre: cost 8 for mismatch
			cost += (((diff >> 4) | (diff >> 5)) & 0x1) << 1 //one out: cost 2 for mismatch
			cost += (((diff >> 1) | diff) & 0x1) << 1
			ds[i] = cost
			start++
		}
	} else if m.k == 6 {
		for i := 0; i < len(ds); i++ {
			diff := kmers[start] ^ a
			cost := (((diff >> 4) | (diff >> 5)) & 0x1) << 2 //centres: cost 4 for mismatch
			cost += (((diff >> 6) | (diff >> 7)) & 0x1) << 2
			cost += (((diff >> 2) | (diff >> 3)) & 0x1) << 1
			cost += (((diff >> 8) | (diff >> 9)) & 0x1) << 1
			cost += (((diff >> 1) | diff) & 0x1)    //edges: cost 1 for mismatch
			cost += (((diff >> 10) | diff>>11) & 0x1)
			ds[i] = cost
			start++
		}
	}
}

type editDistance struct {
	seqs [][]uint16
	rcs  []bool
	k int
	mismatchCost uint16
	insertCost uint16
	deleteCost uint16
}

func NewEditDistance(k int,mismatchCost,insertCost,deleteCost uint16) Measure {
	sm := editDistance{k:k,mismatchCost:mismatchCost,insertCost:insertCost,deleteCost:deleteCost}
	return &sm
}
func (m *editDistance) SetSequences(seqs [][]uint16, rc []bool) {
	m.seqs = seqs
	m.rcs = rc
}
func (m *editDistance) GetSequences() ([][]uint16, []bool) {
	return m.seqs,m.rcs
}
func (m *editDistance) GetSequenceLen(index int) int {
	return len(m.seqs[index])
}
func (m *editDistance) Distances(a uint16, seq, start int, ds []uint16) {
	kmers := m.seqs[seq]
	mismatchCost := m.mismatchCost
	deleteCost := m.deleteCost //delete means the read (in kmers) is missing one from the consensus (a)
	insertCost := m.insertCost
	k := uint16(m.k)

	end := start + len(ds)
	if end > len(kmers) {
		f := len(ds) + len(kmers) - end
		if f < 0 {
			f = 0
		}
		for i := f; i < len(ds); i++ {
			ds[i] = k*mismatchCost
		}
		ds = ds[:f]
	}

	for i := 0; i < len(ds); i++ {
		nextK := kmers[start]
		if nextK == a {
			ds[i] = 0
			start++
			continue
		}
		diff := nextK ^ a
		//collapse the diffs to the odd bits
		diff = diff | (diff >> 1)
		//count diffs in both directions
		var dLHS uint16
		var dRHS uint16
		for j := uint16(0); j < k && ((diff >> (j*2)) & 0x1) == 0; j++ {
			dRHS++
		}
		//if 1 error in diff, return the value
		if dRHS >= k-1 {
			ds[i] = mismatchCost
			start++
			//fmt.Println("1 error (",dRHS," matches)")
			continue
		}
		for j := int(k-1); j >= 0 && ((diff >> uint(j*2)) & 0x1) == 0; j-- {
			dLHS++
		}
		if dLHS+dRHS >= k-1 {
			ds[i] = mismatchCost
			start++
			//fmt.Println("1 error (",dRHS+dLHS," matches)",sequence.KmerString(int(a),int(k)),sequence.KmerString(int(nextK),int(k)))
			continue
		}

		minCost := (k - (dLHS+dRHS))*mismatchCost

		//otherwise, look at diff + right diff (delete found by pushing things right)
		// then right diff + diff (an insert found by pulling things in on the left)
		// and left diff + diff (delete found by pushing things left)
		// and diff + left diff (an insert found by pulling things in on the right)
		rightDiff := ((nextK >> 2) ^ a)
		rightDiff = rightDiff | (rightDiff >> 1)
		leftDiff := ((nextK << 2) ^ a) >> 2
		leftDiff = leftDiff | (leftDiff >> 1)
		//test two deletes first, where nextK was pulled apart
		var rRHS uint16
		var lLHS uint16
		for j := uint16(0); j < k-1 && ((rightDiff >> (j*2)) & 0x1) == 0; j++ { //note: one less to check
			rRHS++
		}
		for j := int(k-2); j >= 0 && ((leftDiff >> uint(j*2)) & 0x1) == 0; j-- {
			lLHS++
		}
		if (dLHS+rRHS >= k-1 || lLHS+dRHS >= k-1) && deleteCost < minCost {
			ds[i] = deleteCost
			start++
			continue
		}
		cost := (k-(dLHS+rRHS))*deleteCost
		if cost < minCost {
			minCost = cost
		}
		cost = (k-(lLHS+dRHS))*deleteCost
		if cost < minCost {
			minCost = cost
		}

		//then the inserts
		var rLHS uint16
		var lRHS uint16
		for j := int(k-2); j >= 0 && ((rightDiff >> uint(j*2)) & 0x1) == 0; j-- {
			rLHS++
		}
		for j := uint16(0); j < k-1 && ((leftDiff >> (j*2)) & 0x1) == 0; j++ {
			lRHS++
		}
		if (dLHS+lRHS >= k-1 || rLHS+dRHS >= k-1) && insertCost < minCost {
			ds[i] = insertCost
			start++
			continue
		}
		cost = (k-(rLHS+dRHS))*insertCost
		if cost < minCost {
			minCost = cost
		}
		cost = (k-(dLHS+lRHS))*insertCost
		if cost < minCost {
			minCost = cost
		}

		//finally, count the mismatches in diff, use that if lower cost
		var misMatches uint16
		for j := uint16(0); j < k; j++ {
			misMatches += (diff >> (j*2)) & 0x1
		}
		if misMatches*mismatchCost < minCost {
			ds[i] = misMatches*mismatchCost
		} else {
			ds[i] = minCost
		}
		start++
	}
}

type matrixDistance struct {
	seqs [][]uint16
	rcs  []bool
	k int
	matrix [][]uint8
}

func NewMatrixDistance(k int,matrix [][]uint8) Measure {
	sm := matrixDistance{k:k,matrix:matrix}
	return &sm
}
func (m *matrixDistance) SetSequences(seqs [][]uint16, rc []bool) {
	m.seqs = seqs
	m.rcs = rc
}
func (m *matrixDistance) GetSequences() ([][]uint16, []bool) {
	return m.seqs,m.rcs
}
func (m *matrixDistance) GetSequenceLen(index int) int {
	return len(m.seqs[index])
}
func (m *matrixDistance) Distances(a uint16, seq, start int, ds []uint16) {
	kmers := m.seqs[seq]
	end := start + len(ds)
	if end > len(kmers) {
		f := len(ds) + len(kmers) - end
		if f < 0 {
			f = 0
		}
		for i := f; i < len(ds); i++ {
			ds[i] = 15
		}
		ds = ds[:f]
	}

	for i := 0; i < len(ds); i++ {
		nextK := kmers[start]
		ds[i] = uint16(m.matrix[a][nextK])
		start++
	}
}
