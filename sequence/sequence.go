package sequence

import (
	"fmt"
)

type Sequence interface {
	GetID() int
	setID(int)
	GetName() string
	String() string
	SubSequence(int, int) Sequence
	StitchSequence(Sequence, int, int) int
	ReverseComplement() Sequence
	Len() int
	GetOffset() int
	GetInset() int
	Kmers(int) []uint64
	ShortKmers(int) []uint16
	Quality() []byte
	Detach() //ensure this sequence does not have slice references on shared arrays
	getBytes() []byte
}

//byteSequence uses an internal mapping of N=0, A=1, B=2, C=3, T=4 storing
//the sequence as a byte slice
type byteSequence struct {
	data    []byte
	quality []byte
	id      int
	offset  int //if a subsequence of a larger sequence
	inset   int
	name    *string
}

//packedSequence stores 4 bases per byte, trading CPU overhead for memory efficiency
type packedSequence struct {
	data     []byte
	quality  []byte
	finalLen int
}

func NewByteSequence(id int, seq string, name string) Sequence {
	data := make([]byte, len(seq), len(seq))
	for i, s := range seq {
		b := byte(s)
		data[i] = ((b >> 1) ^ ((b & 4) >> 2)) & 3
	}
	s := byteSequence{data: data, quality: nil, id: id, offset: 0, inset: 0, name: &name}
	return &s
}

func NewByteSequenceFromKmers(id int, kmers []uint16, k int) Sequence {
	data := make([]byte, len(kmers)+k-1, len(kmers)+k-1)
	for i := 0; i < k-1; i++ {
		data[i] = byte((kmers[0] >> (2 * uint(k-i-1))) & 3)
	}
	for i, kmer := range kmers {
		data[i+k-1] = byte(kmer & 3)
	}
	s := byteSequence{data: data, quality: nil, id: id, offset: 0, inset: 0}
	return &s
}

func (s *byteSequence) StitchSequence(rhs Sequence, minOverlap, maxOverlap int) int {
	if maxOverlap > s.Len() {
		maxOverlap = s.Len()
	}
	if maxOverlap > rhs.Len() {
		maxOverlap = rhs.Len()
	}
	//highest absolute matches
	overlap := 0
	count := 0.0
	as := s.data[s.Len()-maxOverlap:]
	bs := rhs.getBytes()[:maxOverlap]
	for i := minOverlap; i < maxOverlap; i++ {
		c := 0
		for j := 0; j < i; j++ {
			if as[maxOverlap-j-1] == bs[j] {
				c++
			}
		}
		r := float64(c) / float64(i)
		if r > count {
			count = r
			overlap = i
		}
	}
	//and stitch it
	s.data = append(s.data, rhs.getBytes()[overlap:]...)
	//add quality scores if available
	if s.quality != nil {
		q := rhs.Quality()
		if q != nil {
			s.quality = append(s.quality, q[overlap:]...)
		}
	}
	return overlap
}

func (s *byteSequence) ReverseComplement() Sequence {
	bs := make([]byte, len(s.data), len(s.data))
	for i, b := range s.data {
		bs[len(bs)-1-i] = b ^ 3
	}
	//TODO: reverse the quality
	rc := byteSequence{data: bs, quality: nil, id: s.id, offset: s.inset, inset: s.offset, name: s.name}
	return &rc
}

func (s *byteSequence) GetID() int {
	return s.id
}

func (s *byteSequence) setID(id int) {
	s.id = id
}

func (s *byteSequence) GetName() string {
	if s.name == nil {
		return fmt.Sprint(s.id)
	}
	return *(s.name)
}

func (s *byteSequence) String() string {
	buf := make([]byte, len(s.data), len(s.data))
	for i, b := range s.data {
		if b == 0 {
			buf[i] = byte('A')
		} else if b == 1 {
			buf[i] = byte('C')
		} else if b == 2 {
			buf[i] = byte('G')
		} else {
			buf[i] = byte('T')
		}
	}
	return string(buf)
}

func (s *byteSequence) SubSequence(start, end int) Sequence {
	if end > len(s.data) {
		end = len(s.data)
	}
	ss := byteSequence{data: s.data[start:end], id: s.id, offset: s.offset + start, inset: s.inset + len(s.data) - end, name: s.name}
	if s.quality != nil {
		ss.quality = s.quality[start:end]
	}
	return &ss
}

func (s *byteSequence) Detach() {
	newData := make([]byte, len(s.data), len(s.data))
	copy(newData, s.data)
	s.data = newData
	if s.quality != nil {
		newQuality := make([]byte, len(s.quality), len(s.quality))
		copy(newData, s.quality)
		s.quality = newQuality
	}
}

func (s *byteSequence) Len() int {
	return len(s.data)
}

func (s *byteSequence) GetOffset() int {
	return s.offset
}
func (s *byteSequence) GetInset() int {
	return s.inset
}

func (s *byteSequence) Kmers(k int) []uint64 {
	length := len(s.data) - k + 1
	kmers := make([]uint64, length, length)

	var mask uint64
	var v uint64
	for i := 0; i < k; i++ {
		mask = (mask << 2) | 3
		v = (v << 2) | uint64(s.data[i])
	}
	for i := k; i < len(s.data); i++ {
		//store and shift
		kmers[i-k] = v
		v = ((v << 2) | uint64(s.data[i])) & mask
	}
	kmers[length-1] = v
	return kmers
}

//ShortKmers produces a list of kmers (up to 8-mers in length), collapsing homopolymers
func (s *byteSequence) ShortKmers(k int) []uint16 {
	length := len(s.data) - k + 1
	kmers := make([]uint16, length, length)

	var mask uint16
	var v uint16
	for i := 0; i < k; i++ {
		mask = (mask << 2) | 3
		v = (v << 2) | uint16(s.data[i])
	}
	var prev uint16
	index := 0
	for i := k; i < len(s.data); i++ {
		if v != prev || index == 0 {
			//store and shift
			kmers[index] = v
			prev = v
			index++
		}
		v = ((v << 2) | uint16(s.data[i])) & mask
	}
	kmers[index] = v
	index++
	return kmers[:index]
}

func (s *byteSequence) Quality() []byte {
	return s.quality
}

func (s *byteSequence) getBytes() []byte {
	return s.data
}

func KmerString(value int, k int) string {
	bs := make([]byte, k, k)
	for i := k - 1; i >= 0; i-- {
		base := value & 3
		switch base {
		case 0:
			bs[i] = 'A'
		case 1:
			bs[i] = 'C'
		case 2:
			bs[i] = 'G'
		case 3:
			bs[i] = 'T'
		}
		value = value >> 2
	}
	return string(bs)
}
