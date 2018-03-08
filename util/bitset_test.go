package util

import (
	"testing"
)

func Test1CountIntersection(test *testing.T) {
	setA := NewIntSet()
	setB := NewIntSet()
	var count uint
	for i := 1001; i < 3000; i += 5 {
		setA.Add(uint(i))
	}
	for j := 101; j < 2013; j += 3 {
		setB.Add(uint(j))
		if setA.Contains(uint(j)) {
			count++
		}
	}
	if setA.CountIntersection(setB) != count || setB.CountIntersection(setA) != count {
		test.Error("Bad intersection counts:",setA.CountIntersection(setB),setB.CountIntersection(setA),"should be",count)
	}
	next := setB.CountIntersectionTo(setA,int(count+10))
	if next != count {
		test.Error("Bad first intersection to +10 counts (to",count+10,"):",next,"should be",count)
	}
	next = setA.CountIntersectionTo(setB,int(count+10))
	if next != count {
		test.Error("Bad second intersection to +10 counts (to",count+10,"):",next,"should be",count)
	}
	/*if setA.CountIntersectionTo(setB,int(count)) != count || setB.CountIntersectionTo(setA,int(count)) != count {
		test.Error("Bad intersection to count:",setA.CountIntersectionTo(setB,int(count)),setB.CountIntersectionTo(setA,int(count)),"should be",count)
	}
	if setA.CountIntersectionTo(setB,int(count)-2) < count-2 || setB.CountIntersectionTo(setA,int(count)-2) < count-2 {
		test.Error("Bad intersection to count-2:",setA.CountIntersectionTo(setB,int(count)-2),setB.CountIntersectionTo(setA,int(count)-2),"should be",count-2)
	}*/
}
func Test2SharedIDs(test *testing.T) {
	sets := make([]*IntSet, 20)
	for i, _ := range sets {
		sets[i] = NewIntSet()
	}
	counts := make([]int, 500)
	c16 := 0
	c8 := 0
	c4 := 0
	c2 := 0
	for i := 0; i < len(counts); i++ {
		if i % 7 == 0 {
			counts[i] = 16
			c16++
		} else if i % 5 == 0 {
			counts[i] = 8
			c8++
		} else if i % 3 == 0 {
			counts[i] = 4
			c4++
		} else if i % 2 == 0 {
			counts[i] = 2
			c2++
		}
		for j := 0; j < counts[i]; j++ {
			sets[j].Add(uint(i))
		}
	}
	a := GetSharedIDs(sets, 16)
	if len(a) != c16 {
		test.Error("Incorrect number of ids found:",len(a),"vs extected",c16)
	}
	for _, id := range a {
		if counts[id] != 16 {
			test.Error("Incorrect id",id," in",counts[id],"sets but should be in 16")
		}
	}
	a = GetSharedIDsFast(sets, 16)
	if len(a) != c16 {
		test.Error("Incorrect (fast) number of ids found:",len(a),"vs extected",c16)
		test.Error(a)
	}
	for _, id := range a {
		if counts[id] != 16 {
			test.Error("Incorrect (fast) id",id," in",counts[id],"sets but should be in 16")
		}
	}
	a = GetSharedIDs(sets, 8)
	if len(a) != c8+c16 {
		test.Error("Incorrect number of ids found:",len(a),"vs extected",c8)
	}
	for _, id := range a {
		if counts[id] < 8 {
			test.Error("Incorrect id",id," in",counts[id],"sets but should be in 8")
		}
	}
	a = GetSharedIDsFast(sets, 8)
	if len(a) != c8+c16 {
		test.Error("Incorrect (fast) number of ids found:",len(a),"vs extected",c8)
		test.Error(a)
	}
	for _, id := range a {
		if counts[id] < 8 {
			test.Error("Incorrect (fast) id",id," in",counts[id],"sets but should be in 8")
		}
	}
}
