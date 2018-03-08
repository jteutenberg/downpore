package sequence

import (
	"testing"
)

func seqString() string {
	return "GGGAAGTGACTGCCTTAAAATGAGGGTTACCCCTTTTAGTTGACAAGACGCTTGCGGCTATTATGGCTAG"
}

func kmerSet(s string, k int) ([]bool, int) {
	length := 1
	for i := 0; i < k; i++ {
		length = length << 2
	}
	count := 0
	//every 5th and 6th
	ks := make([]bool, length)
	for i := 4; i < len(s)-k-1; i += 5 {
		ks[KmerValue(s[i:i+k])] = true
		ks[KmerValue(s[i+1:i+1+k])] = true
		count += 2
	}
	x := KmerValue(s[0:k])
	if !ks[x] {
		ks[x] = true
		count++
	}
	x = KmerValue(s[1 : k+1])
	if !ks[x] {
		ks[x] = true
		count++
	}
	x = KmerValue(s[len(s)-k:])
	if !ks[x] {
		ks[x] = true
		count++
	}
	return ks, count
}

func Test1Lengths(test *testing.T) {
	s := seqString()
	s1 := NewByteSequence(0, s, nil)
	s2 := NewPackedSequence(0, s, nil)
	if s1.Len() != s2.Len() {
		test.Error("Mismatching lengths: ", s1.Len(), s2.Len(), "from input", len(s))
	}
	for i := 1; i < 5; i++ {
		s1 = NewByteSequence(0, s[:len(s)-i], nil)
		s2 = NewPackedSequence(0, s[:len(s)-i], nil)
		if s1.Len() != s2.Len() {
			test.Error("Mismatching lengths: ", s1.Len(), s2.Len(), "from input", len(s)-i)
		}
	}
}

func Test2String(test *testing.T) {
	s := seqString()
	s1 := NewByteSequence(0, s, nil).String()
	s2 := NewPackedSequence(0, s, nil).String()
	if s1 != s2 {
		test.Error("Mismatching output: \n", s, "\n", s1, "\n", s2)
	}
	for i := 1; i < 5; i++ {
		s1 = NewByteSequence(0, s[:len(s)-i], nil).String()
		s2 = NewPackedSequence(0, s[:len(s)-i], nil).String()
		if s1 != s2 {
			test.Error("Mismatching output: \n", s[:len(s)-i], "\n", s1, "\n", s2,"length",len(s)-i)
		}
	}
}

func Test3Complement(test *testing.T) {
	s := seqString()
	s1 := NewByteSequence(0, s, nil).ReverseComplement().String()
	s2 := NewPackedSequence(0, s, nil).ReverseComplement().String()
	if s1 != s2 {
		test.Error("Mismatching output: \n", s, "\n", s1, "\n", s2)
	}
}
func Test4Subs(test *testing.T) {
	s := seqString()
	s1 := NewByteSequence(0, s, nil)
	s2 := NewPackedSequence(0, s, nil)
	for i := 15; i < 20; i++ {
		a1 := s1.SubSequence(i-15, i).String()
		b1 := s2.SubSequence(i-15, i).String()
		if a1 != b1 {
			test.Error("Mismatching output: \n", s[i-15:i], "\n", a1, "\n", b1)
		}
		a2 := s1.SubSequence(i, i+30).String()
		b2 := s2.SubSequence(i, i+30).String()
		if a2 != b2 {
			test.Error("Mismatching output: \n", s[i:i+30], "\n", a2, "\n", b2)
		}
	}
}
func Test5KmerAt(test *testing.T) {
	s := seqString()
	s1 := NewByteSequence(0, s, nil)
	s2 := NewPackedSequence(0, s, nil)
	for i := 0; i < len(s)-6; i++ {
		a := s1.KmerAt(i, 6)
		b := s2.KmerAt(i, 6)
		if a != b {
			test.Error("Mismatching kmer at: ", i, a, b, s[i:i+6], KmerString(a, 6), KmerString(b, 6))
		}
	}
}

func Test6CountKmers(test *testing.T) {
	s := seqString()
	ks, count := kmerSet(s, 6)
	mask := 0
	for i := 0; i < 6; i++ {
		mask = (mask << 2) | 3
	}
	s1 := NewByteSequence(0, s, nil)
	s2 := NewPackedSequence(0, s, nil)
	c1 := s1.CountKmers(100,6, mask, ks)
	c2 := s2.CountKmers(100,6, mask, ks)
	if c1 != c2 || c1 != count {
		test.Error("Mismatching kmer counts : ", count, c1, c2)
	}
	c1 = s1.CountKmers(7,6, mask, ks)
	c2 = s2.CountKmers(7,6, mask, ks)
	if c1 < 7 || c2 < 7 {
		test.Error("Mismatching partial kmer counts (7):", c1, c2)
	}

	mask = 0
	for i := 0; i < 8; i++ {
		mask = (mask << 2) | 3
	}
	s1 = s1.SubSequence(7, s1.Len()-7)
	s2 = s2.SubSequence(7, s2.Len()-7)
	ks, count = kmerSet(s[7:len(s)-7], 8)
	c1 = s1.CountKmers(100,8, mask, ks)
	c2 = s2.CountKmers(100,8, mask, ks)
	if c1 != c2 || c1 != count {
		test.Error("Mismatching kmer counts : ", count, c1, c2, "\n", s1.String(), "\n", s2.String())
	}
}
func Test7IterateKmers(test *testing.T) {
	s := seqString()
	mask := 0
	for i := 0; i < 6; i++ {
		mask = (mask << 2) | 3
	}
	s1 := NewByteSequence(0, s, nil)
	s2 := NewPackedSequence(0, s, nil)
	for i := 0; i < len(s)-6; i++ {
		k1 := s1.KmerAt(i, 6)
		k2 := s2.KmerAt(i, 6)
		if k1 != k2 {
			test.Error("Mismatching kmers extracted at: ", i, k1, k2)
		}
		k1 = s1.NextKmer(k1, mask, i+6)
		k2 = s2.NextKmer(k2, mask, i+6)
		if k1 != k2 {
			test.Error("Mismatching kmers iterating at: ", i, k1, k2, KmerString(k1, 6), KmerString(k2, 6), " : ", s[i+1:i+7])
		}
	}
}
func Test8Segments(test *testing.T) {
	s := seqString()
	mask := 0
	for i := 0; i < 6; i++ {
		mask = (mask << 2) | 3
	}
	s1 := NewByteSequence(0, s, nil)
	s2 := NewPackedSequence(0, s, nil)
	ks, count := kmerSet(s, 6)
	seg1 := make([]int, count*2+2)
	seg2 := make([]int, count*2+2)
	s1.WriteSegments(seg1, 6, mask, ks)
	s2.WriteSegments(seg2, 6, mask, ks)
	err := false
	for i, x := range seg1 {
		err = err || (x != seg2[i])
	}
	if err {
		test.Error("Mismatching seed seq: \n", seg1, "\n", seg2)
	}
	for i, _ := range seg1 {
		seg1[i] = 4321
		seg2[i] = 4321
	}
	s1.SubSequence(2, s1.Len()-2).WriteSegments(seg1, 6, mask, ks)
	s2.SubSequence(2, s2.Len()-2).WriteSegments(seg2, 6, mask, ks)
	err = false
	for i, x := range seg1 {
		err = err || (x != seg2[i])
	}
	if err {
		test.Error("Mismatching seed seq: \n", seg1, "\n", seg2)
	}
}

func Test9Packing(test *testing.T) {
	s := "CGGT"
	data := make([]byte, 2) //expect an empty second byte
	packBytes([]byte(s), data)
	if data[0] != 0x6B {
		test.Error("Incorrect byte encoding:", data[0], "for", s)
	}
	if data[1] != 0 {
		test.Error("Overflow:", data)
	}

	s = "CGGTCGGTCGGTCGGTCGGT"
	data = make([]byte, len(s)/4+1)
	packBytes([]byte(s), data)
	for i, b := range data[:len(data)-1] {
		if b != 0x6B {
			test.Error("Incorrect byte encoding at:", i, data[i], "for", s)
		}
	}
	if data[len(data)-1] != 0 {
		test.Error("Overflow:", data)
	}
}

func Test10ShortKmers(test *testing.T) {
	s := seqString()
	s1 := NewByteSequence(0, s, nil)
	s2 := NewPackedSequence(0, s, nil)

	k1 := s1.ShortKmers(6, false)
	k2 := s2.ShortKmers(6, false)
	if len(k1) != len(k2) {
		test.Error("Incorrect kmer seq lengths:", len(k1), len(k2), "\n", k1, "\n", k2)
	}
	err := false
	for i, x := range k1 {
		err = err || x != k2[i]
	}
	if err {
		test.Error("Incorrect kmer seqs:\n", k1, "\n", k2)
	}
	k1 = s1.ShortKmers(3, true)
	k2 = s2.ShortKmers(3, true)
	if len(k1) != len(k2) {
		test.Error("Incorrect kmer seq lengths:", len(k1), len(k2), "from sequences length", s1.Len(), s2.Len(), "\n", k1, "\n", k2)
	}
	err = false
	for i, x := range k1 {
		err = err || x != k2[i]
	}
	if err {
		test.Error("Incorrect kmer seqs:\n", k1, "\n", k2)
	}
}
