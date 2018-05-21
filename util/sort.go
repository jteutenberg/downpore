package util

import (
	"sort"
)

type ByValue struct {
	Values []int
	IDs []int
}
func (d *ByValue) Len() int {
	return len(d.IDs)
}
func (d *ByValue) Less(i, j int) bool {
	return d.Values[i] < d.Values[j]
}
func (d *ByValue) Swap(i, j int) {
	d.IDs[i], d.IDs[j] = d.IDs[j], d.IDs[i]
	d.Values[i], d.Values[j] = d.Values[j], d.Values[i]
}

func SortByValue(ids, values []int) {
	sort.Sort(&ByValue{IDs:ids,Values:values})
}
