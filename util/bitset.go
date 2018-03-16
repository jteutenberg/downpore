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
	set := IntSet{make([]uint64, 50), 1, 0, 0}
	return &set
}

func NewIntSetCapacity(capacity int) *IntSet {
	set := IntSet{make([]uint64, capacity/64+1), 1, 0, 0}
	return &set
}

func NewIntSetFromInts(values []int) *IntSet {
	var max int
	for _, v := range values {
		if v > max {
			max = v
		}
	}
	set := IntSet{make([]uint64, max/64+1), 1, 0, 0}
	s := &set
	s.AddInts(values)
	return s
}

func NewIntSetFromUInts(values []uint) *IntSet {
	var max uint
	for _, v := range values {
		if v > max {
			max = v
		}
	}

	set := IntSet{make([]uint64, max/64+1), 1, 0, 0}
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
	if int(index) >= len(set.vs) {
		newVs := make([]uint64, index+1)
		copy(newVs, set.vs)
		set.vs = newVs
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

func (set *IntSet) CountIntersectionTo(other *IntSet, maxCount int) uint {
	start := set.start
	end := set.end
	if other.start > start {
		start = other.start
	}
	if end > other.end {
		end = other.end
	}

	/*count := 0
	for ; start <= end && count < maxCount; start++ {
		count += bits.OnesCount64(set.vs[start] & other.vs[start])
	}
	return uint(count)*/
	return countIntersectionToAsm(set.vs[start:end+1], other.vs[start:end+1], maxCount)
}

func countIntersectionToAsm(a, b []uint64, maxCount int) uint

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
		newVs := make([]uint64, set.end+1)
		copy(newVs, set.vs)
		set.vs = newVs
	}
	for i := start; i <= end; i++ {
		set.vs[i] |= other.vs[i]
	}
}

//gets union, and first three soft-unions (up to 4+ shared sets)
func getSoftUnion4Asm(vs []uint64) (uint64, uint64, uint64, uint64)

//gets soft unions with 5-8+ shared sets
func getSoftUnion8Asm(vs []uint64) (uint64, uint64, uint64, uint64)

//gets soft union with 16+ shared sets
func getSoftUnion16Asm(vs []uint64) (uint64, uint64, uint64, uint64)

func getSoftUnion8(vs []uint64) (uint64, int) {
	zeroCount := 0 //empty vs
	n := len(vs)
	maxZeroes := n - 8
	var v, v2, v3, v4, v5, v6, v7, v8 uint64
	for j := 0; j < n && zeroCount <= maxZeroes; j++ {
		m := vs[j]
		v8 |= v7 & m
		v7 |= v6 & m
		v6 |= v5 & m
		v5 |= v4 & m
		v4 |= v3 & m
		v3 |= v2 & m
		v2 |= v & m
		v |= m
		if m == 0 {
			zeroCount++
		}
	}
	return v8, zeroCount
}

func getSoftUnion16(vs []uint64) (uint64, int) {
	zeroCount := 0 //empty vs
	n := len(vs)
	maxZeroes := n - 16
	var v, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16 uint64
	for j := 0; j < n && zeroCount <= maxZeroes; j++ {
		m := vs[j]
		v16 |= v15 & m
		v15 |= v14 & m
		v14 |= v13 & m
		v13 |= v12 & m
		v12 |= v11 & m
		v11 |= v10 & m
		v10 |= v9 & m
		v9 |= v8 & m
		v8 |= v7 & m
		v7 |= v6 & m
		v6 |= v5 & m
		v5 |= v4 & m
		v4 |= v3 & m
		v3 |= v2 & m
		v2 |= v & m
		v |= m
		if m == 0 {
			zeroCount++
		}
	}
	return v16, zeroCount
}

func GetSharedIDs(sets []*IntSet, minCount int, fast bool) (ids []uint) {
	if minCount > 24 {
		fast = false
	}
	n := len(sets)
	start := uint(len(sets[0].vs))
	end := uint(0)
	lens := make([]uint, n, n)
	vs := make([][]uint64, n)
	for i := 0; i < n; i++ {
		vs[i] = sets[i].vs
		lens[i] = uint(sets[i].end + 1)
		if sets[i].start < start {
			start = sets[i].start
		}
		if sets[i].end > end {
			end = sets[i].end
		}
	}
	nextVs := make([]uint64, n)
	for i := start; i <= end; i++ {
		for j := 0; j < n; j++ {
			if lens[j] > i {
				nextVs[j] = vs[j][i]
			} else {
				nextVs[j] = 0
			}
		}
		var v uint64 //union
		if minCount >= 13 {
			v13, v14, v15, v16 := getSoftUnion16Asm(nextVs)
			if minCount >= 16 {
				v = v16
			} else if minCount == 15 {
				v = v15
			} else if minCount == 14 {
				v = v14
			} else {
				v = v13
			}
		} else if minCount >= 8 {
			v5, v6, v7, v8 := getSoftUnion8Asm(nextVs)
			if minCount >= 8 {
				v = v8
			} else if minCount == 7 {
				v = v7
			} else if minCount == 6 {
				v = v6
			} else {
				v = v5
			}
		} else {
			v1, v2, v3, v4 := getSoftUnion4Asm(nextVs)
			if minCount == 4 {
				v = v4
			} else if minCount == 3 {
				v = v3
			} else if minCount == 2 {
				v = v2
			} else {
				v = v1
			}
		}
		//extract ids from v
		if v != 0 {
			if ids == nil {
				ids = make([]uint, 0, 20)
			}
			if fast {
				var shifted uint
				for v != 0 {
					zs := uint(bits.TrailingZeros64(v))
					ids = append(ids, (i<<6)+shifted+zs)
					v = v >> (zs + 1)
					shifted += zs + 1
				}
			} else {
				ids = addSoftUnionIDs(v, nextVs, minCount, ids, uint(i<<6))
			}
		}
	}
	return ids
}

/*
//GetSharedIDs IDs that are contained by at least
//minCount of the provided sets
func GetSharedIDs(sets []*IntSet, minCount int) []uint {
	if minCount > 16 && minCount < 24 {
		minCount = 16
	}
	start := uint(len(sets[0].vs))
	end := uint(0)
	n := len(sets)
	lens := make([]uint, n, n)
	vs := make([][]uint64, n)
	var ids []uint
	for i := 0; i < n; i++ {
		vs[i] = sets[i].vs
		lens[i] = uint(sets[i].end+1) //uint(len(sets[i].vs))
		if sets[i].start < start {
			start = sets[i].start
		}
		if sets[i].end > end {
			end = sets[i].end
		}
	}
	maxZeroes := n - minCount
	for i := start; i <= end; i++ {
		count := 0
		var v uint64  //union
		if minCount >= 16 {
			v, count = getSoftUnion16(vs,lens,i)
		} else if minCount >= 8 {
			v, count = getSoftUnion8(vs,lens,i)
			minCount = 8
		} else {
			var v2 uint64 //all bits that appear 2+ times
			var v3 uint64 // "" 3+ times
			var v4 uint64 // "" 4+ times
			//find how many zeroed longs we have here, build up the union as we go
			for j := 0; j < n && count <= maxZeroes; j++ {
				if lens[j] > i {
					n := vs[j][i]
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
		}
		//TODO: report bits directly from the appropriate pseudo-union?
		if count <= maxZeroes && v != 0 { //enough non-zero longs
			//we now have the union of the sets. Check each ID
			bit := Bit
			zs := uint(bits.TrailingZeros64(v))
			bit <<= zs
			v >>= zs
			for j := zs; j < 64 && v != 0; j++ {
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
}*/

//validate each set bit in v against the original vs
func addSoftUnionIDs(v uint64, vs []uint64, minCount int, ids []uint, offset uint) []uint {
	bit := Bit
	zs := uint(bits.TrailingZeros64(v))
	bit <<= zs
	v >>= zs
	n := len(vs)
	for j := zs; j < 64 && v != 0; j++ {
		if (Bit & v) != 0 { //a member of the union
			count := 0
			//find which sets contained this one (make sure we ignore shortened vs slices)
			for k := 0; k < n; k++ {
				if (vs[k] & bit) != 0 {
					count++
					if count >= minCount { //this bit has appeared in enough sets
						if ids == nil {
							ids = make([]uint, 0, 20)
						}
						ids = append(ids, uint(offset+j))
						break
					}
				} else if n-k+count <= minCount {
					break
				}
			}
		}
		v >>= 1
		bit <<= 1
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
	return true, (index << 6) + subIndex + zs
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
