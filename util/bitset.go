package util

import (
	"fmt"
	"math/bits"
)

const (
	Bit     uint64 = 1
	AllBits uint64 = 0xFFFFFFFFFFFFFFFF
)

type IntSet struct {
	vs    []uint64
	start uint
	end   uint
	count uint //size is no greater than this
}

func NewIntSet() *IntSet {
	set := IntSet{make([]uint64, 0, 50), 1, 0, 0}
	return &set
}

func NewIntSetFromInts(values []int) *IntSet {
	set := IntSet{make([]uint64, 0, 50), 1, 0, 0}
	s := &set
	s.AddInts(values)
	return s
}

func NewIntSetFromUInts(values []uint) *IntSet {
	set := IntSet{make([]uint64, 0, 50), 1, 0, 0}
	s := &set
	for _, v := range values {
		s.Add(v)
	}
	return s
}

func (set *IntSet) AddInts(values []int) {
	for _, v := range values {
		set.Add(uint(v))
	}
}

func (set *IntSet) Contains(x uint) bool {
	index := x >> 6
	if index < set.start || index > set.end {
		return false
	}
	subIndex := x & 0x3F
	return (set.vs[index] & (Bit << subIndex)) != 0
}

func (set *IntSet) Add(x uint) {
	index := x >> 6
	subIndex := x & 0x3F
	bit := Bit << subIndex
	for int(index) >= len(set.vs) {
		set.vs = append(set.vs, 0)
	}
	old := set.vs[index]
	if (old & bit) != 0 { //already exists
		return
	}
	set.vs[index] = old | bit
	if index < set.start || set.start > set.end {
		set.start = index
	}
	if index > set.end {
		set.end = index
	}
	set.count++
}

func (set *IntSet) Remove(x uint) {
	index := x >> 6
	if index > set.end {
		return
	}
	subIndex := x & 0x3F
	bit := Bit << subIndex
	old := set.vs[index]
	if (old & bit) == 0 {
		return //nothing to remove
	}
	set.vs[index] ^= bit
	if index == set.start || index == set.end {
		set.reduceStartEnd()
	}
	set.count--
}

func (set *IntSet) reduceStartEnd() {
	for set.start <= set.end && set.vs[set.start] == 0 {
		set.start++
	}
	for set.end >= set.start && set.vs[set.end] == 0 {
		set.end--
	}
	if set.start > set.end {
		set.start = uint(len(set.vs)) + 1
		set.end = 0
	}
}

func (set *IntSet) IsEmpty() bool {
	return set.start > set.end
}

func (set *IntSet) Clear() {
	for set.start <= set.end {
		set.vs[set.start] = 0
		set.start++
	}
	set.end = 0
	set.start = uint(len(set.vs)) + 1
	set.count = 0
}

func (set *IntSet) GetFirstID() (bool, uint) {
	if set.IsEmpty() {
		return false, 0
	}
	v := set.vs[set.start]
	return true, set.start*64 + uint(bits.TrailingZeros64(v))
}

func (set *IntSet) CountIntersection(other *IntSet) uint {
	start := set.start
	end := set.end
	if other.start > start {
		start = other.start
	}
	if end > other.end {
		end = other.end
	}
	count := 0
	for ; start <= end; start++ {
		count += bits.OnesCount64(set.vs[start] & other.vs[start])
	}
	return uint(count)
}

func (set *IntSet) Intersect(other *IntSet) {
	for ; set.start < other.start && set.start <= set.end; set.start++ {
		set.vs[set.start] = 0
	}
	for ; set.end > other.end && set.end >= set.start; set.end-- {
		set.vs[set.end] = 0
	}
	for i := set.start; i <= set.end; i++ {
		set.vs[i] &= other.vs[i]
	}
	set.reduceStartEnd()
}

func (set *IntSet) RemoveAll(other *IntSet) {
	start := set.start
	end := set.end
	if other.start > start {
		start = other.start
	}
	if end > other.end {
		end = other.end
	}
	for i := start; i <= end; i++ {
		set.vs[i] &= (^other.vs[i])
	}
	set.reduceStartEnd()
}

func (set *IntSet) Union(other *IntSet) {
	start := other.start
	end := other.end
	if start < set.start {
		set.start = start
	}
	if end > set.end {
		set.end = end
		for int(set.end) >= len(set.vs) {
			set.vs = append(set.vs, 0)
		}
	}
	for i := start; i <= end; i++ {
		set.vs[i] |= other.vs[i]
	}
}

//GetSharedIDs IDs that are contained by at least
//minCount of the provided sets
func GetSharedIDs(sets []*IntSet, minCount int) []uint {
	start := uint(len(sets[0].vs))
	end := uint(0)
	n := len(sets)
	lens := make([]uint, n, n)
	var ids []uint
	for i := 0; i < n; i++ {
		lens[i] = uint(len(sets[i].vs))
		if sets[i].start < start {
			start = sets[i].start
		}
		if sets[i].end > end {
			end = sets[i].end
		}
	}
	maxZeroes := n - minCount
	for i := start; i <= end; i++ {
		var v uint64  //union
		var v2 uint64 //all bits that appear 2+ times
		var v3 uint64 // "" 3+ times
		var v4 uint64 // "" 4+ times
		//find how many zeroed longs we have here, build up the union as we go
		count := 0
		for j := 0; j < n && count <= maxZeroes; j++ {
			if lens[j] > i {
				n := sets[j].vs[i]
				v4 |= v3 & n
				v3 |= v2 & n
				v2 |= v & n
				v |= n
				if n == 0 {
					count++
				}
			} else {
				count++
			}
		}
		if minCount >= 3 {
			if minCount >= 4 {
				v = v4
			} else {
				v = v3
			}
		} else {
			v = v2
		}
		//TODO: report bits directly from the appropriate pseudo-union?
		if count <= maxZeroes && v != 0 { //enough non-zero longs
			//we now have the union of the sets. Check each ID
			bit := Bit
			zs := uint(bits.TrailingZeros64(v))
			bit <<= zs
			v >>= zs
			for j := zs; /*uint(0);*/ j < 64 && v != 0; j++ {
				if (Bit & v) != 0 { //a member of the union
					count = 0
					//find which sets contained this one (make sure we ignore shortened vs slices)
					for k := 0; k < n; k++ {
						if lens[k] > i && (sets[k].vs[i]&bit) != 0 {
							count++
							if count >= minCount { //this bit has appeared in enough sets
								if ids == nil {
									ids = make([]uint, 0, 20)
								}
								ids = append(ids, uint(i*64+j))
								break
							}
						} else if k-count > maxZeroes {
							break
						}
					}
				}
				v >>= 1
				bit <<= 1
			}
		}
	}
	return ids
}

func (set *IntSet) GetNextID(x uint) (bool, uint) {
	x++
	index := x >> 6
	if index > set.end {
		return false, 0
	}
	subIndex := uint(0)
	if index < set.start {
		index = set.start
	} else {
		subIndex = x & 0x3F
		filter := AllBits << subIndex
		if (filter & set.vs[index]) == 0 {
			index++
			subIndex = 0
		}
	}
	for index <= set.end && set.vs[index] == 0 {
		index++
	}
	if index > set.end {
		return false, 0
	}
	v := set.vs[index] >> subIndex
	zs := uint(bits.TrailingZeros64(v))
	return true, index*64 + subIndex + zs
}

func (set *IntSet) AsInts() []int {
	//consider inlining the loop to speed this up
	ids := make([]int, 0, set.count)
	for ok, id := set.GetFirstID(); ok; ok, id = set.GetNextID(id) {
		ids = append(ids, int(id))
	}
	return ids
}

func (set *IntSet) CountMembers() uint {
	count := 0
	for i := set.start; i <= set.end; i++ {
		count += bits.OnesCount64(set.vs[i])
	}
	set.count = uint(count)
	return set.count
}

//Size gets an upper bound on the size of this set. After a call to CountMembers this value
//is accurate until elements are removed.
func (set *IntSet) Size() uint {
	return set.count
}

func (set *IntSet) String() string {
	s := "{"
	first := true
	for ok, v := set.GetFirstID(); ok; ok, v = set.GetNextID(v) {
		if first {
			first = false
			s = fmt.Sprint(s, v)
		} else {
			s = fmt.Sprint(s, ",", v)
		}
	}
	return s + "}"
}
