package alignment

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
