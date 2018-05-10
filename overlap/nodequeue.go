//Largely taken from the documentation for the heap package

// This example demonstrates a priority queue built using the heap interface.
package overlap

import (
	"container/heap"
)

//QueueItem is something we manage in a priority queue.
type QueueItem struct {
	value    *Node
	rc	bool //relative reverse-complementedness for this queue
	distance int //lower is wanted earlier
	index int // The index of the item in the heap.
}

// A PriorityQueue implements heap.Interface and holds Items.
type PriorityQueue []*QueueItem

func NewPriorityQueue() PriorityQueue {
	states := make(PriorityQueue, 0, 20)
	heap.Init(&states)
	return states
}

func (pq PriorityQueue) Len() int { return len(pq) }
func (pq PriorityQueue) Less(i, j int) bool {
	return pq[i].distance < pq[j].distance
}

func (pq PriorityQueue) Swap(i, j int) {
	pq[i], pq[j] = pq[j], pq[i]
	pq[i].index = i
	pq[j].index = j
}

func (pq *PriorityQueue) PushNode(n *Node, distance int, rc bool) {
	item := QueueItem{value:n, distance:distance, rc:rc}
	heap.Push(pq, &item)
}

func (pq *PriorityQueue) PopNode() (*Node,int,bool) {
	item := heap.Pop(pq).(*QueueItem)
	return item.value,item.distance,item.rc
}

func (pq *PriorityQueue) Push(x interface{}) {
	n := len(*pq)
	item := x.(*QueueItem)
	item.index = n
	*pq = append(*pq, item)
}

func (pq *PriorityQueue) Pop() interface{} {
	old := *pq
	n := len(old)
	item := old[n-1]
	item.index = -1 // for safety
	*pq = old[0 : n-1]
	return item
}

