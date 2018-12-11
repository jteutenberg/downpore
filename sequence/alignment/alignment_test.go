package alignment

import (
	"fmt"
	"github.com/jteutenberg/downpore/sequence"
	"testing"
	//"math/rand"
)

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

func Test1EditDistance(test *testing.T) {
	m := NewEditDistance(5)
	truth := []string{"AACAA","CCCAA","CCAAA","CACAC","CAAAA","AAAAC","ACCCA","ACCAA","AACCA"}
	truthmer := [][]uint16{make([]uint16,len(truth))}
	for i, t := range truth {
		truthmer[0][i] = uint16(sequence.KmerValue(t))
	}
	m.SetSequences(truthmer,nil)
	inserters := []string{"","T","TT","TTT"}

	ds := make([]uint16,1)
	for j, t := range truth {
		//2x6 inserts (shifting left and right)
		for insert_size := 1; insert_size < 3; insert_size++ {
			for i := 0; i <= len(t); i++ {
				var insert string
				if i <= len(t)-insert_size {
					insert = t[:i]+ inserters[insert_size] + t[i:len(t)-insert_size] //LHS unchanged
					m.Distances(uint16(sequence.KmerValue(insert)), 0, j, ds)
					if (insert_size == 1 && ds[0] != 1) || (insert_size > 1 && ds[0] <= 1) {
						test.Error("Edit distance ",ds[0]," for ",t," vs ",insert)
					}
				}
				if i > insert_size {
					insert = t[insert_size:i] + inserters[insert_size] + t[i:] //LHS shifted
					m.Distances(uint16(sequence.KmerValue(insert)), 0, j, ds)
					if (insert_size == 1 && ds[0] != 1) || (insert_size > 1 && ds[0] <= 1) {
						test.Error("Edit distance ",ds[0]," for ",t," vs ",insert)
					}
				}
			}
		}
		//2x5 deletes (adding a T at each end)
		for insert_size := 1; insert_size < 3; insert_size++ {
			for i := 0; i <= len(t)-insert_size; i++ {
				deleted := inserters[insert_size] + t[:i] + t[i+insert_size:]
				m.Distances(uint16(sequence.KmerValue(deleted)), 0, j, ds)
				if (insert_size == 1 && ds[0] != 1) || (insert_size > 1 && ds[0] <= 1) {
					test.Error("Edit distance ",ds[0]," for ",t," vs ",deleted)
				}
				deleted = t[:i] + t[i+insert_size:]+inserters[insert_size]
				m.Distances(uint16(sequence.KmerValue(deleted)), 0, j, ds)
				if (insert_size == 1 && ds[0] != 1) || (insert_size > 1 && ds[0] <= 1) {
					test.Error("Edit distance ",ds[0]," for ",t," vs ",deleted)
				}
			}
		}
		//5 mismatches
	}
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
