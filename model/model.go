package model

import (
	"bufio"
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/sequence"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

type Model interface {
	GetK() uint
	SetSequences([][]uint16, []bool)
	GetSequences() ([][]uint16, []bool)
	GetSequenceLen(int) int
	//Distance is the relative to the expected current difference between two kmers
	Distance(a, b uint16) uint16
	Distances(a uint16, sequence, start int, ds []uint16)
	//DistanceRC is the difference between the reverse complements of a and b
	DistanceRC(a, b uint16) uint16
	//Distance2D is the sum of the differences between a and b and the difference between the reverse-complements
	Distance2D(a, b uint16) uint16
	Clone() Model
}

type model struct {
	k uint
	is2D bool
	levels []uint16
	rcLevels []uint16
	seqs [][]uint16 //sequences as squiggles
	rcSeqs [][]uint16 //sequences as reverse complement squiggles
	originals [][]uint16
	rc []bool
}

func NewModel(filename string, is2D bool) Model {
	m := model{0,is2D,nil,nil, nil, nil, nil, nil}
	fin, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer fin.Close()

	a := byte('A')
	t := byte('T')
	bin := bufio.NewReader(fin)
	var levels []float64
	for buf, berr := bin.ReadBytes('\n'); len(buf) > 0 || berr == nil; buf, berr = bin.ReadBytes('\n') {
		if buf[0] < a || buf[0] > t {
			continue
		}
		//split into two tokens
		s := string(buf[:len(buf)-1])
		tokens := strings.Split(s, "\t")
		seq := sequence.NewByteSequence(0,tokens[0],nil)
		v := uint16(seq.ShortKmers(seq.Len(),true)[0])
		//the first sequence.. sort out storage
		if m.k == 0 {
			m.k = uint(seq.Len())
			count := 1 << (2*m.k)
			levels = make([]float64, count, count)
		}
		levels[v], _ = strconv.ParseFloat(tokens[1], 64)
	}
	m.levels = make([]uint16, len(levels), len(levels))
	m.rcLevels = make([]uint16, len(levels), len(levels))
	//sort the levels, make the 20th to 80th span 100 values. Threshold.
	temp := make([]float64, len(levels), len(levels))
	copy(temp, levels)
	sort.Float64s(temp)
	minLevel := temp[len(temp)/5]
	maxLevel := temp[len(temp) - len(temp)/5]
	maxDelta := maxLevel-minLevel //threshold anything above this
	f := 255.0/maxDelta
	for i, level := range levels {
		scaled := (level-temp[0])*f
		if scaled > 10000.0 {
			scaled = 10000.0
		}
		m.levels[i] = uint16(scaled)
		m.rcLevels[seeds.ReverseComplement(uint(i),m.k)] = m.levels[i]
	}
	return &m
}

func (m *model) Clone() Model {
	nm := model{k:m.k, levels:m.levels, rcLevels:m.rcLevels}
	return &nm
}
func (m *model) GetK() uint {
	return m.k
}

func (m *model) GetSequences() ([][]uint16, []bool) {
	return m.originals, m.rc
}

func (m *model) GetSequenceLen(index int) int {
	return len(m.originals[index])
}

func (m *model) SetSequences(seqs [][]uint16, rc []bool) {
	m.originals = seqs
	m.rc = rc
	m.seqs = make([][]uint16, len(seqs), len(seqs))
	m.rcSeqs = make([][]uint16, len(seqs), len(seqs))
	for i, seq := range seqs {
		s := make([]uint16, len(seq), len(seq))
		rs := make([]uint16, len(seq), len(seq))
		for j, b := range seq {
			s[j] = m.levels[b]
			rs[j] = m.rcLevels[b]
		}
		m.seqs[i] = s
		m.rcSeqs[i] = rs
	}
}

func (m *model) Distances(a uint16, sequence, start int, ds []uint16) {
	if m.is2D {
		m.distances2D(a, sequence, start, ds)
		return
	}
	var level uint16
	var levels []uint16
	if m.rc[sequence] {
		level = m.rcLevels[a]
		levels = m.rcSeqs[sequence]
	} else {
		level = m.levels[a]
		levels = m.seqs[sequence]
	}
	kmers := m.originals[sequence]
	end := start+len(ds)
	if end > len(kmers) {
		for i := len(ds)+ len(kmers)-end; i < len(ds); i++ {
			ds[i] = 1000
		}
		ds = ds[:len(ds) - end + len(kmers)]
		end = len(kmers)
	}
	var d uint16
	var b uint16
	for i := 0; i < len(ds); i++ {
		b = levels[start]
		if b < level {
			d = 1+level-b
		} else if b > level {
			d = 1+b - level
		} else if a == kmers[start] {
			d = 0
		} else {
			d = 1
		}
		if d > 50 {
			ds[i] = 50
		} else {
			ds[i] = d
		}
		start++
	}
}
func (m *model) distances2D(a uint16, sequence, start int, ds []uint16) {
	level := m.levels[a]
	rcLevel := m.rcLevels[a]
	kmers := m.originals[sequence]
	levels := m.seqs[sequence]
	rcLevels := m.rcSeqs[sequence]
	end := start+len(ds)
	if end > len(kmers) {
		for i := len(ds)+ len(kmers)-end; i < len(ds); i++ {
			ds[i] = 1000
		}
		ds = ds[:len(ds) - end + len(kmers)]
		end = len(kmers)
	}
	var d uint16
	var b uint16
	var rcb uint16
	for i := 0; i < len(ds); i++ {
		if a == kmers[start] {
			ds[i] = 0
			continue
		}
		b = levels[start]
		rcb = rcLevels[start]
		if b < level {
			d = 1+level-b
		} else if b > level {
			d = 1+b - level
		} else {
			d = 1
		}
		if rcb < rcLevel {
			d += 1+rcLevel-rcb
		} else if rcb > rcLevel {
			d += 1+rcb - rcLevel
		} else {
			d += 1
		}
		d /= 2
		if d > 50 {
			d = 50
		}
		ds[i] = d
		start++
	}
}

func (m *model) Distance(a, b uint16) uint16 {
	if a == b {
		return 0
	}
	a = m.levels[a]
	b = m.levels[b]
	var d uint16
	if a < b {
		d = b-a
	} else {
		d = a-b
	}
	if d >= 49 {
		return 50
	}
	return 1+d
}
func (m *model) DistanceRC(a, b uint16) uint16 {
	if a == b {
		return 0
	}
	a = m.rcLevels[a]
	b = m.rcLevels[b]
	var d uint16
	if a < b {
		d = b-a
	} else {
		d = a-b
	}
	if d >= 49 {
		return 50
	}
	return 1+d
}
func (m *model) Distance2D(a, b uint16) uint16 {
	return m.Distance(a,b) + m.DistanceRC(a,b)
}
