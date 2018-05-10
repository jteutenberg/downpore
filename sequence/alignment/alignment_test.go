package alignment

import (
	"fmt"
	//"github.com/jteutenberg/downpore/sequence"
	"testing"
	//"math/rand"
)

/*
type simpleMeasure struct {}

func (m *simpleMeasure) Distance(a, b uint16) uint16 {
	if a == b {
		return 0
	}
	return 3
}
func (m *simpleMeasure) Distances(a uint16, sequence, start int, ds []uint16) {
	end := start+len(ds)
		ds[i] = m.Distance(a,b)
	}
}

type fivemerMeasure struct{}
func (m *fivemerMeasure) Distance(a, b uint16) uint16 {
	if a == b {
		return 0
	}
	if ((a ^ b) & 252) == 0 {
		return 2
	}
	return 3
}
func (m *fivemerMeasure) DistanceRC(a, b uint16) uint16 {
	return 0
}
func (m *fivemerMeasure) Distance2D(a, b uint16) uint16 {
	return 0
}
func (m *fivemerMeasure) DistancesRC(a uint16, bs []uint16, ds []uint16) {
	return
}
func (m *fivemerMeasure) Distances(a uint16, bs, ds []uint16) {
	for i,b := range bs {
		ds[i] = m.Distance(a,b)
	}
}

func makeDTW(numSeq int, k int) Aligner {
	var m Measure
	if k == 5 {
		m = &(fivemerMeasure{})
	} else {
		m = &(simpleMeasure{})
	}
	return NewDTWAligner(10,2, m, 15, k)
}

func errorise(seq string, k, add, remove, perm int, random *rand.Rand) []uint16 {
	bases := []string{"A","C","G","T"}
	ns := ""
	for i, s := range seq {
		r := random.Intn(100)
		if i < 2 {
			r = 100
		}
		if r < add+perm {
			ns = ns+bases[random.Intn(4)]
		} else if r >= add+perm+remove {
			ns = ns+string(s)
		}
	}
	return sequence.NewByteSequence(0,ns).ShortKmers(k)
}

func checkOutput(dtw Aligner, ss [][]uint16, k int, ans []uint16, t *testing.T) {
	rcs := make([]bool, len(ss), len(ss))
	kmers,costs, _ := dtw.GlobalConsensus(ss,rcs)
	cs := ""
	i := 0
	bad := false
	for kmer := range kmers {
		cost := <-costs
		fmt.Print(sequence.KmerString(int(kmer), k),",")
		cs = fmt.Sprint(cs,",",cost)
		if i < len(ans) && kmer != ans[i] {
			t.Error("Mismatch at index",i,"answer",sequence.KmerString(int(ans[i]),k),"!=",sequence.KmerString(int(kmer),k))
			bad = true
		}
		i++
	}
	if i < len(ans)-2 { //allow 2 at the end in case of modified sequences
		t.Error("Incomplete call:",i,"kmers instead of",len(ans))
	}
	if bad {
		t.Error("\nCorrect sequence:\n",sequence.NewByteSequenceFromKmers(-1,ans,k).String())
		t.Error("\nInput:")
		for _,s := range ss {
			t.Error("\n",sequence.NewByteSequenceFromKmers(-1,s,k).String())
		}
	}
	fmt.Println()
	fmt.Println(cs)
}
*/
func Test0Asm(test *testing.T) {
	//before := []uint16{65535,0,65535,65535,65535,65535,65535,65535,65535,0,65535,65535,65535,2,2,5,0,1,65535,65535,65535,65535,65535,65535,65535,65535,65535,65535,65535,0,65535,65535}
	//before := []uint16{32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,0,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767}
	//dist := []uint16{1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,3,0,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3}
	before := []uint16{32767,32767,0,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767,32767}
	dist := []uint16{50,50,5,50,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256,1256}
	/*dist := make([]uint16, 32,32)
	for i, _ := range dist {
		if i %2 == 0 {
			dist[i] = 8
		} else {
			dist[i] = 0
		}
	}*/
	output := make([]uint16, 32, 32)
	cost := updateOffsetsAsm(dist, before, output, 300)
	for _, x := range before {
		if x > 100 {
			fmt.Print("- ")
		} else {
			fmt.Print(x," ")
		}
	}
	fmt.Println()
	fmt.Println(dist)
	for _, x := range output{
		if x > 100 {
			fmt.Print("- ")
		} else {
			fmt.Print(x," ")
		}
	}
	fmt.Println("cost :",cost)
}
/*
func Test1SingleSeq(test *testing.T) {
	dtw := makeDTW(1,3)
	s := sequence.NewByteSequence(0, "AAACAAGTGCGTGAGAAGCCTTT")
	ks := s.ShortKmers(3)
	ss := make([][]uint16, 1, 1)
	ss[0] = ks
	checkOutput(dtw,ss,3,ks,test)
}

func Test2SameSeqs(test *testing.T) {
	dtw := makeDTW(1,3)
	s := sequence.NewByteSequence(0, "AAACAAGTGCGTGAGAAGCCTTT")
	ks := s.ShortKmers(3)
	ss := make([][]uint16, 3, 3)
	ss[0], ss[1], ss[2] = ks,ks,ks
	checkOutput(dtw,ss,3,ks,test)
}
func Test3SimpleSeqs(test *testing.T) {
	dtw := makeDTW(3,3)
	k1 := sequence.NewByteSequence(0, "AAATAAGTGCGTGAGTAGCCTTT").ShortKmers(3)
	k2 := sequence.NewByteSequence(0, "AAACAAGAGCGTGAGAAGCTTTT").ShortKmers(3)
	k3 := sequence.NewByteSequence(0, "AAACAAGTGCATGAGAAGCCTTA").ShortKmers(3)
	ss := make([][]uint16, 3, 3)
	ss[0], ss[1], ss[2] = k1,k2,k3
	checkOutput(dtw,ss,3,sequence.NewByteSequence(0, "AAACAAGTGCGTGAGAAGCCTTT").ShortKmers(3),test)
}

func Test4Inserts(test *testing.T) {
	dtw := makeDTW(3,3)
	k1 := sequence.NewByteSequence(0, "TTAATACAAGTGCGTGAGAAGCACTTT").ShortKmers(3)
	k2 := sequence.NewByteSequence(0, "TTAAACAAGTGGCGTGAGAAGCCTTT").ShortKmers(3)
	k3 := sequence.NewByteSequence(0, "TTAAACAAGTGCGTGACGAAGCCTTT").ShortKmers(3)
	ss := make([][]uint16, 3, 3)
	ss[0], ss[1], ss[2] = k1,k2,k3
	checkOutput(dtw,ss,3, sequence.NewByteSequence(0, "TTAAACAAGTGCGTGAGAAGCCTTT").ShortKmers(3), test)
}
func Test4Deletes(test *testing.T) {
	dtw := makeDTW(3,3)
	k1 := sequence.NewByteSequence(0, "TTAACAAGTGCGTGAGAAGCTTT").ShortKmers(3)
	k2 := sequence.NewByteSequence(0, "TTAAACAAGGCGTGAGAAGCCTTT").ShortKmers(3)
	k3 := sequence.NewByteSequence(0, "TTAAACAAGTGCGTAGAAGCCTTT").ShortKmers(3)
	ss := make([][]uint16, 3, 3)
	ss[0], ss[1], ss[2] = k1,k2,k3
	checkOutput(dtw,ss,3, sequence.NewByteSequence(0, "TTAAACAAGTGCGTGAGAAGCCTTT").ShortKmers(3), test)
}

func Test5Edges(test *testing.T) {
	dtw := makeDTW(5,3)
	k1 := sequence.NewByteSequence(0, "AAACAAGTGCGTGAGAAGCCTTT").ShortKmers(3)
	k2 := sequence.NewByteSequence(0,      "AGTGCGTGAGAAGCCT"  ).ShortKmers(3)
	k3 := sequence.NewByteSequence(0,   "ACAAGTGCGTGAGAAGC"    ).ShortKmers(3)
	ss := make([][]uint16, 3, 3)
	ss[0], ss[1], ss[2] = k1,k2,k3
	//checkOutput(dtw,ss,3,sequence.NewByteSequence(0, "AGTGCGTGAGAAGC").ShortKmers(3), test)
	//checkOutput(dtw,ss,3,sequence.NewByteSequence(0, "ACAAGTGCGTGAGAAGCCT").ShortKmers(3), test)
	//checkOutput(dtw,ss,3,sequence.NewByteSequence(0, "AAACAAGTGCGTGAGAAGCCTTT").ShortKmers(3), test)
	checkOutput(dtw,ss,3,sequence.NewByteSequence(0, "AAACAAGTGCGTGAGAAGC").ShortKmers(3), test)
}
func Test6Deep(test *testing.T) {
	depth := 10
	dtw := makeDTW(depth,5)
	r := rand.New(rand.NewSource(551255))
	seq := "TAATAAGTGCGTGAGTAGCCTTTAAATAAGTGCGTGAGTAGCCTTT"
	ss := make([][]uint16, depth, depth)
	for i := 0; i < len(ss); i++ {
		ss[i] = errorise(seq,5, 6, 6, 6, r)
	}
	checkOutput(dtw,ss,5,sequence.NewByteSequence(0,seq).ShortKmers(5),test)
}
func Test7Lengths(test *testing.T) {
	dtw := makeDTW(3,3)
	seq := "TAATAAGTGCGTGAGTAGCCTTTAAATAAGTGCGTGAGTAGCCTTAA"
	seq = seq + seq + seq + seq
	ss := make([][]uint16, 4, 4)
	ss[0] = sequence.NewByteSequence(0, seq).ShortKmers(3)
	ss[1] = sequence.NewByteSequence(0, seq).ShortKmers(3)
	step := 6
	ns := seq[:step/2]
	ps := ""
	i := 0
	for ; i < len(seq)-step/2-step; i+=step {
		ps = ps + seq[i:i+step]+"G"
		ns = ns + seq[i+step/2:i+step/2+step-1]
	}
	ps = ps + seq[i:]
	ns = ns + seq[i:]
	ss[2] = sequence.NewByteSequence(0,ps).ShortKmers(3)
	ss[3] = sequence.NewByteSequence(0,ns).ShortKmers(3)
	checkOutput(dtw,ss,3,sequence.NewByteSequence(0,seq).ShortKmers(3),test)
}
*/
