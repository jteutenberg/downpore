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
	ReverseComplement() Sequence
	Len() int
	GetOffset() int
	GetInset() int
	KmerAt(int, int) int
	NextKmer(int, int, int) int
	ShortKmers(int, bool) []uint16
	Quality() []byte
	SetQuality([]byte)
	Detach() //ensure this sequence does not have slice references on shared arrays
	CountKmers(int, int, int, []bool) int
	CountKmersBetween(int, int,int, int, int, []bool) int
	WriteSegments([]int, int, int, []bool)
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

//packedSequence stores 4 bases per byte
type packedSequence struct {
	data     []byte
	quality  []byte
	id       int
	offset   int
	inset    int
	name     *string
	length   int //length of sequence in bases
	firstLen int //number of bases in first byte
	finalLen int //number of bases in final byte
}

func NewByteSequence(id int, seq string, name *string) Sequence {
	data := make([]byte, len(seq), len(seq))
	for i, s := range seq {
		b := byte(s)
		data[i] = ((b >> 1) ^ ((b & 4) >> 2)) & 3
	}
	s := byteSequence{data: data, quality: nil, id: id, offset: 0, inset: 0, name: name}
	return &s
}

func packBytes(seq []byte, data []byte)

func NewPackedSequence(id int, seq string, name *string) Sequence {
	length := len(seq) / 4
	internalLength := length * 4
	finalLength := len(seq) - internalLength //partial bytes at the end
	data := make([]byte, (len(seq)+3)/4)
	if internalLength >= 4 {
		packBytes([]byte(seq[:internalLength]), data)
	}
	// and the partial one at the end (with trailing zeroes)
	if finalLength > 0 {
		var b byte
		for i := finalLength; i > 0; i-- {
			nb := byte(seq[len(seq)-i])
			nb = ((nb >> 1) ^ ((nb & 4) >> 2)) & 3
			b = (b << 2) | nb
		}
		if finalLength < 4 {
			b = b << uint(8-finalLength*2)
		}
		data[len(data)-1] = b
	}
	s := packedSequence{data: data, quality: nil, id: id, offset: 0, inset: 0, name: name, length: len(seq), firstLen: 4, finalLen: finalLength}
	if s.finalLen > s.length {
		s.finalLen = s.length
	}
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

//Find determines the best matching position of a subsequence within start to end of this sequence
/*func (s *byteSequence) Find(subs Sequence, subsKmers *util.IntSet, k int) int {
	//build kmers of this sequence
	kmers := util.NewIntSet()
	ks := s.ShortKmers(k)
	for _, n := range ks {
		subKmers.Add(uint(n))
	}
	//find those shared
	kmers.Intersect(subsKmers)
	//generate gapped-seed sequences
	//and get the best linear match

}*/

func (s *byteSequence) ReverseComplement() Sequence {
	bs := make([]byte, len(s.data), len(s.data))
	for i, b := range s.data {
		bs[len(bs)-1-i] = b ^ 3
	}
	var qs []byte
	if s.quality != nil {
		qs = make([]byte, len(s.quality), len(s.quality))
		for i, q := range s.quality {
			qs[len(qs)-1-i] = q
		}
	}
	rc := byteSequence{data: bs, quality: qs, id: s.id, offset: s.inset, inset: s.offset, name: s.name}
	return &rc
}

func (s *packedSequence) ReverseComplement() Sequence {
	//gross reverse
	bs := make([]byte, len(s.data), len(s.data))
	for i, b := range s.data {
		//byte complement
		b = ^b
		//and byte reverse
		bs[len(bs)-1-i] = ((b & 3) << 6) | ((b & 12) << 2) | ((b & 48) >> 2) | ((b & 192) >> 6)
	}
	var qs []byte
	if s.quality != nil {
		qs = make([]byte, len(s.quality), len(s.quality))
		for i, q := range s.quality {
			// reverse for quality
			qs[len(qs)-1-i] = q
		}
	}
	rc := packedSequence{data: bs, quality: qs, id: s.id, offset: s.inset, inset: s.offset, firstLen: s.finalLen, finalLen: s.firstLen, name: s.name, length: s.length}
	return &rc
}

func (s *byteSequence) GetID() int {
	return s.id
}
func (s *packedSequence) GetID() int {
	return s.id
}

func (s *byteSequence) setID(id int) {
	s.id = id
}
func (s *packedSequence) setID(id int) {
	s.id = id
}

func (s *byteSequence) GetName() string {
	if s.name == nil {
		return fmt.Sprint(s.id)
	}
	return *(s.name)
}
func (s *packedSequence) GetName() string {
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
func (s *packedSequence) String() string {
	buf := make([]byte, len(s.data)*4, len(s.data)*4)
	j := s.firstLen*2 - 2
	var count uint
	for _, b := range s.data[:len(s.data)-1] {
		for j >= 0 {
			buf[count] = (b >> uint(j)) & 3
			count++
			j -= 2
		}
		j = 6
	}
	b := s.data[len(s.data)-1]
	last := 8 - s.finalLen*2
	for j >= last {
		buf[count] = (b >> uint(j)) & 3
		count++
		j -= 2
	}
	for i, b := range buf {
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
	return string(buf[:s.length])
}

func (s *byteSequence) CountKmers(upTo, k, mask int, kmers []bool) int {
	seed := s.KmerAt(0, k) >> 2
	count := 0
	for i := k - 1; i < len(s.data); i++ {
		seed = s.NextKmer(seed, mask, i)
		if kmers[seed] {
			count++
			if count >= upTo {
				break
			}
		}
	}
	return count
}

func (s *byteSequence) CountKmersBetween(from, to, upTo, k, mask int, kmers []bool) int {
	seed := s.KmerAt(from, k) >> 2
	count := 0
	for i := from+k - 1; i < to; i++ {
		seed = s.NextKmer(seed, mask, i)
		if kmers[seed] {
			count++
			if count >= upTo {
				break
			}
		}
	}
	return count
}

func (s *byteSequence) WriteSegments(segments []int, k, mask int, seeds []bool) {
	seed := s.KmerAt(0, k) >> 2
	kmerIndex := 0
	prev := 0
	count := 0
	for i := k - 1; i < len(s.data); i++ {
		seed = s.NextKmer(seed, mask, i)
		if seeds[seed] {
			segments[count] = kmerIndex - prev
			segments[count+1] = seed
			prev = kmerIndex + k
			count += 2
		}
		kmerIndex++
	}
	segments[count] = len(s.data) - prev
}

func packedWriteSegments(data []byte, skipFront, skipBack, k int, seeds []bool, segments []int)
func packedCountKmers(data []byte, upTo, skipFront, skipBack, k int, seeds []bool) int

func (s *packedSequence) CountKmers(upTo, k, mask int, seeds []bool) int {
	return packedCountKmers(s.data, upTo, 4-s.firstLen, 4-s.finalLen, k, seeds)
}
func (s *packedSequence) CountKmersBetween(from, to,upTo, k, mask int, seeds []bool) int {
	//shrink in from the left and right a bit, i.e. let us undercount sometimes
	start := (from+4-s.firstLen + 3)/4
	end := (to+4-s.firstLen)/4
	return packedCountKmers(s.data[start:end], upTo, 4-s.firstLen, 4-s.finalLen, k, seeds)
}
func (s *packedSequence) WriteSegments(segments []int, k, mask int, seeds []bool) {
	packedWriteSegments(s.data, 4-s.firstLen, 4-s.finalLen, k, seeds, segments)
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

func (s *packedSequence) SubSequence(start, end int) Sequence {
	if end > s.length {
		end = s.length
	}
	end--                         //make this inclusive
	off := start + 4 - s.firstLen //offset in bases
	offByte := off / 4            //offset index
	off -= offByte * 4            //remainder offset in first byte
	in := end + 4 - s.firstLen    //offset to the final base
	inByte := in / 4
	in -= inByte * 4
	data := s.data[offByte : inByte+1]
	ss := packedSequence{data: data, id: s.id, offset: s.offset + start, inset: s.inset + s.length - end, name: s.name, firstLen: 4 - off, finalLen: in + 1, length: end - start + 1}
	if s.quality != nil {
		ss.quality = s.quality[start : end+1]
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

func (s *packedSequence) Detach() {
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
func (s *packedSequence) Len() int {
	return s.length
}

func (s *byteSequence) GetOffset() int {
	return s.offset
}
func (s *packedSequence) GetOffset() int {
	return s.offset
}

func (s *byteSequence) GetInset() int {
	return s.inset
}
func (s *packedSequence) GetInset() int {
	return s.inset
}

func (s *byteSequence) KmerAt(index, k int) int {
	var v int
	for i := index; i < index+k; i++ {
		v = (v << 2) | int(s.data[i])
	}
	return v
}

func packedKmerAt(data []byte, offset, k int) int32

func (s *packedSequence) KmerAt(index, k int) int {
	return int(packedKmerAt(s.data, index+4-s.firstLen, k))
}

func (s *byteSequence) NextKmer(current, mask int, nextBaseIndex int) int {
	return ((current << 2) | int(s.data[nextBaseIndex])) & mask
}
func (s *packedSequence) NextKmer(current, mask int, nextBaseIndex int) int {
	nextBaseIndex += 4 - s.firstLen
	b := s.data[nextBaseIndex/4]
	subIndex := uint(3-(nextBaseIndex&3)) << 1 //in bits
	b = (b >> subIndex) & 3
	return ((current << 2) | int(b)) & mask
}

//ShortKmers produces a list of kmers (up to 8-mers in length), collapsing homopolymers if required
func (s *byteSequence) ShortKmers(k int, collapse bool) []uint16 {
	length := len(s.data) - k + 1
	kmers := make([]uint16, length)

	var mask uint16
	var v uint16
	for i := 0; i < k; i++ {
		mask = (mask << 2) | 3
		v = (v << 2) | uint16(s.data[i])
	}
	var prev uint16
	index := 0
	for i := k; i < len(s.data); i++ {
		if !collapse || v != prev || index == 0 {
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

func (s *packedSequence) ShortKmers(k int, collapse bool) []uint16 {
	length := s.length - k + 1
	kmers := make([]uint16, length)
	v := s.KmerAt(0, k)
	var mask int
	for i := 0; i < k; i++ {
		mask = (mask << 2) | 3
	}
	var prev int
	index := 0
	for i := k; i < s.length; i++ {
		if !collapse || v != prev || index == 0 {
			//store and shift
			kmers[index] = uint16(v)
			prev = v
			index++
		}
		v = s.NextKmer(v, mask, i)
	}
	kmers[index] = uint16(v)
	index++
	return kmers[:index]
}

func (s *byteSequence) Quality() []byte {
	return s.quality
}
func (s *packedSequence) Quality() []byte {
	return s.quality
}

func (s *byteSequence) SetQuality(q []byte) {
	s.quality = q
}
func (s *packedSequence) SetQuality(q []byte) {
	s.quality = q
}

func KmerValue(s string) int {
	value := 0
	for _, c := range s {
		b := byte(c)
		d := ((b >> 1) ^ ((b & 4) >> 2)) & 3
		value = (value << 2) | int(d)
	}
	return value
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
