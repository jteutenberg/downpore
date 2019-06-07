package overlap

import (
	"fmt"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/util"
	"sync"
)

//Node represents one long sequence and all its overlapping sequences
//Arcs out belong to all sequences that overhang the ends
type Node struct {
	id int
	colour int
	internalSequences []int
	sequences []*sequenceArc
	Consensus sequence.Sequence
	in []*Arc
	out []*Arc
}

//sequenceArc is the connection between overlap contigs (Nodes) and sequences (SequenceNodes)
type sequenceArc struct {
	sequence *SequenceNode
	node *Node
	approximate bool //whether the offset+length information is approximate or accurate
	offset int //offset of the overlap within the sequence
	length int //length of the overlap within the sequence
	rc bool //whether this node is reverse-complement of the sequence
}

//An ordered list of nodes connected by a shared sequence
type SequenceNode struct {
	id int //the sequence's id
	colour int
	length int

	isRC int //votes
	isNotRC int

	nodes []*sequenceArc

	covered bool //whether nodes exist at both ends of this sequence
	coveredFront bool //a node exists at the front
	coveredBack bool //a node exists at the back
}

//Arc is link between two overlaps based on one or more shared sequences
type Arc struct {
	from *Node //the end of this node..
	to *Node // .. overlaps the beginning of this node
	length int //estimated distance in bases between the overlaps
	fromRC bool
	toRC bool
}

type OverlapGraph struct {
	nodes []*Node
	sequences []*SequenceNode //ordered lists of nodes by sequence
	lock *sync.Mutex
	nextColour int
}

func NewOverlapGraph(maxSeqs int) *OverlapGraph {
	var lock sync.Mutex
	g := OverlapGraph{nodes:make([]*Node, 0, 10000), sequences:make([]*SequenceNode, maxSeqs,maxSeqs), lock:&lock, nextColour:2}
	return &g
}

func (n *Node) ArcLength(index int) int {
	return n.out[index].length
}

func (n *Node) removeInArc(a *Arc) {
	i := 0
	for i < len(n.in) && n.in[i] != a {
		i++
	}
	//in[i] is now a
	for i < len(n.in)-1 {
		n.in[i] = n.in[i+1]
		i++
	}
	n.in = n.in[:len(n.in)-1]
}

func (n *Node) removeOutArc(a *Arc) {
	i := 0
	for i < len(n.out) && n.out[i] != a {
		i++
	}
	//out[i] is now a
	for i < len(n.out)-1 {
		n.out[i] = n.out[i+1]
		i++
	}
	n.out = n.out[:len(n.out)-1]
}

//addArcBetween connects two nodes a and b that share at least one sequence.
//A number of temporary IntSets need to be provided (for re-use). These need not be empty.
//nextColour is a forward-colour to apply to unvisited nodes. Sequences are not coloured by this function.
func (g *OverlapGraph) addArcBetween(a *Node, b *Node, aSeq *util.IntSet, bSeq *util.IntSet, aRC *util.IntSet, bRC *util.IntSet, nextColour int) {
	aSeq.Clear()
	bSeq.Clear()
	aRC.Clear()
	bRC.Clear()
	for _, arc := range a.sequences {
		id := uint(arc.sequence.id)
		aSeq.Add(id)
		if arc.rc {
			aRC.Add(id)
		}
	}
	for _, arc := range b.sequences {
		id := uint(arc.sequence.id)
		if aSeq.Contains(id) {
			bSeq.Add(id)
			if arc.rc {
				bRC.Add(id)
			}
		}
	}
	aSeq.Intersect(bSeq)
	aRC.Intersect(bSeq)
	aSeq.RemoveAll(aRC)
	totalShared := bSeq.CountMembers() //at this point bSeq is the set of all shared sequences
	if totalShared == 0 {
		//somehow, these are not actually connected
		fmt.Println("not actually connected..")
		return
	}
	bSeq.RemoveAll(bRC) // and now bSeq is reduced to just the forward sequences in b

	//1. decide whether these nodes are reverse-complement of one another
	same := aRC.CountIntersection(bRC) + aSeq.CountIntersection(bSeq)
	rc := same < totalShared / 2

	//2. if one node already has a colour, apply the appropriate colour to the other node, and reverse it if necessary
	//NOTE: nodes are all assumed to be forward-coloured at this point
	fmt.Println("Adding arc between colours:",a.colour," ",b.colour," RC:",rc)
	if a.colour != 0 && b.colour == 0 {
		if rc {
			b.reverse()
		}
		b.colour = a.colour
	} else if b.colour != 0 && a.colour == 0 {
		if rc {
			a.reverse()
		}
		a.colour = b.colour
	} else if a.colour == 0 && b.colour == 0 {
		//find the majority arc colour in the forward and rc sequences
		c1 := 0 //count of "a is forward"
		c2 := 0 //count of "a is not forward"
		for _, arc := range a.sequences {
			if arc.sequence.colour == 0 {
				continue
			}
			id := uint(arc.sequence.id)
			if aSeq.Contains(id) {
				if arc.sequence.colour == nextColour {
					c1++
				} else {
					c2++
				}
			} else if aRC.Contains(id) {
				if arc.sequence.colour == nextColour {
					c2++
				} else {
					c1++
				}
			}
		}
		//flip the appropriate a or b
		fmt.Println("Counted ",c1,"vs",c2," for a to be RCed")
		if c1 >= c2 {
			//a is forward
			if rc {
				b.reverse()
			}
		} else {
			a.reverse()
			if !rc {
				b.reverse() //to match
			}
		}
		//and assign a colour
		a.colour = nextColour
		b.colour = nextColour
	}

	//3. determine the length and direction of the gap. Add the arc.
	aSeq.Union(aRC) //how holds all shared sequences
	offset := 0
	badCount := 0
	for i := len(a.sequences)-1; i >= 0; i-- {
		arc := a.sequences[i]
		id := uint(arc.sequence.id)
		if aSeq.Contains(id) {
			//find the corresponding arc in b
			for j, barc := range b.sequences {
				bid := uint(barc.sequence.id)
				if bid == id {
					if arc.rc != barc.rc {
						badCount++
						//and remove the arcs (sequences) from both nodes.
						a.sequences[i] = a.sequences[len(a.sequences)-1]
						a.sequences = a.sequences[:len(a.sequences)-1]
						b.sequences[j] = b.sequences[len(b.sequences)-1]
						b.sequences = b.sequences[:len(b.sequences)-1]
						//then the arcs (nodes) from this sequence. These are in order.
						ns := barc.sequence.nodes
						rem := 0 //should add up to 2
						for k := 0; k < len(ns); k++ {
							if ns[k].node == a || ns[k].node == b {
								rem++
							}
							if rem > 0 && rem+k < len(ns) {
								ns[k] = ns[k+rem]
							}
						}
						barc.sequence.nodes = barc.sequence.nodes[:len(barc.sequence.nodes)-rem]
						//fmt.Println("next would be:",barc.offset-arc.offset-arc.length)
						totalShared--
					} else {
						//oldOff := offset
						if arc.rc { // reverse if sequence is RC. barc should be RC too
							offset += barc.offset - arc.offset - arc.length
						} else {
							offset += arc.offset - barc.offset - barc.length
						}
						//fmt.Println("next offset:",offset-oldOff)
					}
					break
				}
			}
		}
	}
	fmt.Println(badCount," arcs want to incorrectly reverse complement a node.",totalShared," do not.")
	if totalShared == 0 {
		return
	}
	offset /= int(totalShared)
	fmt.Println("makes: ",offset)
	if offset < 0 {
		g.addArc(b, a, -offset, false, false)
	} else {
		g.addArc(a, b, offset, false, false)
	}
}

func (g *OverlapGraph) addArc(from *Node, to *Node, size int, fromRC, toRC bool) {
	arc := Arc{from:from, to:to, length:size, fromRC:fromRC, toRC:toRC}
	if to.in == nil {
		to.in = make([]*Arc, 0, 5)
	}
	if from.out == nil {
		from.out = make([]*Arc, 0, 5)
	}
	from.out = append(from.out, &arc)
	to.in = append(to.in, &arc)
	//bubble down so that out arcs are in distance-order
	i := len(from.out)-2
	for ; i >= 0 && size < from.out[i].length; i-- {
		from.out[i+1] = from.out[i]
	}
	from.out[i+1] = &arc
	i = len(to.in)-2
	for ; i >= 0 && size < to.in[i].length; i-- {
		to.in[i+1] = to.in[i]
	}
	to.in[i+1] = &arc
}

//Reverse-complements the consensus sequence and flips the node's colour (if applicable)
//Any sequence arcs in/out have their polarity reversed.
func (n *Node) reverse() {
	n.Consensus = n.Consensus.ReverseComplement()
	if n.colour != 0 {
		n.colour = rcColour(n.colour)
	}
	for _, arc := range n.sequences {
		arc.rc = !arc.rc
	}
}

func rcColour(colour int) int {
	return colour ^ 1
}

func isRCColour(colour int) bool {
	return colour & 1 == 0
}

//isAdjacent checks whether a points to b
func (a *Node) isAdjacent(b *Node) bool {
	for _, arc := range a.out {
		if arc.to == b {
			return true
		}
	}
	return false
}

//Adds a node using the sequence ids, lengths and offsets from the given contig.
func (g *OverlapGraph) AddNode(contig *SeedContig, consensus sequence.Sequence) {
	offsets := contig.Offsets
	g.lock.Lock()
	n := Node{ id:len(g.nodes), colour:0,sequences:make([]*sequenceArc,len(contig.Parts),len(contig.Parts)), Consensus:consensus, in:nil, out:nil}
	g.nodes = append(g.nodes, &n)

	for i, s := range contig.Parts {
		//create an entry for this sequence if necessary
		seq := g.sequences[s]
		if seq == nil {
			seq = &SequenceNode{id:s, length: contig.SeqLengths[i],nodes:make([]*sequenceArc,0,5)}
			g.sequences[s] = seq
		}
		//create the connection between the node and sequence
		arc := sequenceArc{sequence:seq, node:&n, offset:contig.Offsets[i], length:contig.Lengths[i],approximate:contig.Approximate[i],rc:contig.ReverseComplement[i]}
		n.sequences[i] = &arc
		//test sequence coverage
		if arc.offset < arc.length {
			seq.coveredFront = true //close to the beginning of the sequence
			seq.covered = seq.coveredBack
		}
		if arc.offset + arc.length*2 > seq.length {
			seq.coveredBack = true
			seq.covered = seq.coveredFront
		}
		//find where this node should go in the sequence
		index := len(seq.nodes)-1 //index to node before the new one
		for ; index >= 0; index-- {
			if seq.nodes[index].offset < offsets[i] {
				break
			}
		}

		//add this node to the sequence
		seq.nodes = append(seq.nodes, nil) //make room
		//bubble down
		for i := len(seq.nodes)-1; i > index; i-- {
			if i > 0 {
				seq.nodes[i] = seq.nodes[i-1]
			}
		}
		seq.nodes[index+1] = &arc
	}
	g.lock.Unlock()
}

//mergeNodes combines the consensus, takes the union of sequences and updates their offsets/lengths
func (g *OverlapGraph) mergeNodes(a, b *Node) bool {
	//test for sufficient shared sequences
	as := util.NewIntSet()
	bs := util.NewIntSet()
	appa := util.NewIntSet()
	for _, arc := range a.sequences {
		as.Add(uint(arc.sequence.id))
		if !arc.approximate {
			appa.Add(uint(arc.sequence.id))
		}
	}
	for _, arc := range b.sequences {
		if arc.offset > 0 {
			bs.Add(uint(arc.sequence.id))
		}
	}
	as.Intersect(bs)
	if as.CountIntersection(appa) > 1 { //non-trivial overlap
		as.CountMembers()
		fmt.Println("Merging ",a.id,b.id,", sharing ",as.Size()," out of ",len(a.sequences)," and ",len(b.sequences))
		//check for colour merge and possible reverse-complement propagation
		if a.colour != b.colour {
			rcCount := 0
			fwCount := 0
			for _, arc := range a.sequences {
				if as.Contains(uint(arc.sequence.id)) {
					//find the same arc in b, see if they have the same rc-ness
					for _, barc := range b.sequences {
						if barc.sequence.id == arc.sequence.id {
							if arc.rc == barc.rc {
								fwCount++
							} else {
								rcCount++
							}
							break
						}
					}
				}
			}
			if rcCount > fwCount {
				//reverse-complement *every* node with this colour!
				fmt.Println("Need to rc the node.")
			}
		}
		//nodes are now guaranteed to be of the same complementarity

		//estimate the overlap size
		overlapLength := 0
		count := 0
		for _, arc := range a.sequences {
			if !arc.approximate && as.Contains(uint(arc.sequence.id)) {
				for _, barc := range b.sequences {
					if barc.sequence.id == arc.sequence.id {
						if arc.rc {
							fmt.Println(barc.offset+barc.length - arc.offset," C. Overlap in nodes size", arc.length," and ",barc.length, " in sequence ",barc.sequence.id)
							overlapLength += barc.offset+barc.length - arc.offset
						} else {
							fmt.Println(arc.offset+arc.length - barc.offset," overlap in nodes size", arc.length," and ",barc.length, " in sequence ",barc.sequence.id)
							overlapLength += arc.offset+arc.length - barc.offset
						}
						count++
						break
					}
				}
			}
		}
		overlapLength /= count
		fmt.Println("Avg overlap from ",count," = ",overlapLength)
		//combine consensus
		minOverlap := overlapLength-100
		if minOverlap < 0 {
			minOverlap = 5
		}
		
		//overlapLength = a.Consensus.StitchSequence(b.Consensus,minOverlap,overlapLength+100)
		//fmt.Println("Actual overlap is ",overlapLength)
		//fmt.Println("Sizes are: ",a.Consensus.Len(),b.Consensus.Len())
		fmt.Println(a.Consensus.String())
		fmt.Println(b.Consensus.String())
		/*for _, arc := range a.seqs {
			//update each shared sequence
			if as.Contains(uint(arc.sequence.id)) {
			
			}
		}
		for _, arc := range b.seqs {
			//and add new sequences
			if !as.Contains(uint(arc.sequence.id)) {

			}
		}*/
		return true
	}
	return false
}

//AddBridge adds a new overlap that spans a gap between two existing nodes
//The expectation is that this will merge the three nodes into one new contig
func (g *OverlapGraph) AddBridge(contig *SeedContig, consensus sequence.Sequence) {
	//find the before and after nodes
	s := g.sequences[contig.Parts[0]]
	var before *Node
	var after *Node
	for i, n := range s.nodes { //run through the arcs
		if before == nil && n.offset > contig.Offsets[0] {
			fmt.Println("Bridge from prev end:",s.nodes[i-1].offset+s.nodes[i-1].length," bridge start:",contig.Offsets[0]," length ",contig.Lengths[0],". index ",i-1)
			before = s.nodes[i-1].node
		} else if n.offset > contig.Offsets[0]+contig.Lengths[0] {
			after = s.nodes[i-1].node
			fmt.Println(" to bridge end at ",contig.Offsets[0]+contig.Lengths[0]," which should hit next node at ",s.nodes[i-1].offset,". index ",i-1)
			break
		}
	}
	if after == nil {
		fmt.Println("no end found for bridge to ",contig.Offsets[0]+contig.Lengths[0],". Final node at ",len(s.nodes)-1," starts at ",s.nodes[len(s.nodes)-1].offset)
	}
	fmt.Println(before == after)
	g.lock.Lock()
	//create a new node with merged consensus

	//remove two of the old nodes, replace the third

	//add the union of sequences and update their offsets
	
	g.lock.Unlock()
}

func buildContig(seqSet *util.IntSet, leftNode, rightNode *Node, rcBridge bool) *SeedContig {
	contig := SeedContig{Combined:nil,Parts:seqSet.AsInts()}
	contig.ReverseComplement = make([]bool, len(contig.Parts), len(contig.Parts))
	contig.Approximate = make([]bool, len(contig.Parts), len(contig.Parts))
	contig.Offsets = make([]int, len(contig.Parts), len(contig.Parts))
	contig.Lengths = make([]int, len(contig.Parts), len(contig.Parts))
	edgeBuffer := 20 //extra bases to ensure an overlap. 16+ should work.
	for k, id := range contig.Parts {
		i := 0
		j := 0
		for leftNode.sequences[i].sequence.id != id {
			i++
		}
		for rightNode.sequences[j].sequence.id != id {
			j++
		}
		fmt.Println("Adding ",k,"th sequence to the bridge. Offsets:",leftNode.sequences[i].offset,rightNode.sequences[j].offset," with the rc-ness: ",leftNode.sequences[i].rc,rightNode.sequences[j].rc)
		if leftNode.sequences[i].rc != rcBridge { //the bridge will extend from left, so keep its rc-ness
			contig.Offsets[k] = rightNode.sequences[j].offset+rightNode.sequences[j].length - edgeBuffer
			contig.Lengths[k] = leftNode.sequences[i].offset - contig.Offsets[k] + edgeBuffer*2
		} else {
			contig.Offsets[k] = leftNode.sequences[i].offset + leftNode.sequences[i].length - edgeBuffer //from the end of the left..
			contig.Lengths[k] = rightNode.sequences[j].offset - contig.Offsets[k] + edgeBuffer*2
		}
		contig.ReverseComplement[k] = leftNode.sequences[i].rc
		contig.Approximate[k] = false
	}
	return &contig
}


//GetBridgableNodes generates SeedContigs for all segements with sufficient reads and no existing consensus sequence.
func (g *OverlapGraph) GetBridgableContigs(minCoverage uint) []*SeedContig {
	left := util.NewIntSet()
	right := util.NewIntSet()
	bridges := make([]*SeedContig, 0, 100)
	usedBefore := make([]bool, len(g.nodes), len(g.nodes))
	usedAfter := make([]bool, len(g.nodes), len(g.nodes))
	for _, s := range g.sequences {
		if s != nil {
			//look for a gaps between all nodes
			prev := s.nodes[0] //an arc to a node
			for i := 1; i < len(s.nodes); i++ {
				n := s.nodes[i]
				reversed := prev.rc //any old way to collapse the symmetry works: for now, reverse complement on the LHS
				alreadyBridged := (!reversed && (usedAfter[prev.node.id] || usedBefore[n.node.id])) || (reversed && (usedAfter[n.node.id] || usedBefore[prev.node.id]))
				if !alreadyBridged && n.offset > prev.offset+prev.length {
					//test whether or not there are enough sequences spanning the gap
					left.Clear()
					right.Clear()
					for _, a := range prev.node.sequences {
						if !a.approximate {
							left.Add(uint(a.sequence.id))
						}
					}
					for _, a := range n.node.sequences {
						if !a.approximate {
							right.Add(uint(a.sequence.id))
						}
					}
					left.Intersect(right)
					coverage := left.CountMembers()
					//fmt.Println("Testing for a bridge between ",prev.node.id," and ",n.node.id," of distance ",n.offset - prev.offset-prev.length," at offsets ",n.offset,prev.offset,". Coverage:",coverage)
					if coverage > minCoverage {
						contig := buildContig(left, prev.node, n.node, prev.rc)
						bridges = append(bridges, contig)
						if prev.rc { //note: must be the same symmetry removal as above
							usedBefore[prev.node.id] = true
							usedAfter[n.node.id] = true
						} else {
							usedAfter[prev.node.id] = true
							usedBefore[n.node.id] = true
						}
					}
				}
				prev = n
			}
		}
	}
	return bridges
}

func (g *OverlapGraph) GenerateArcs() {
	//check the node->sequence->node lists for important arcs
	for _, seq := range g.sequences {
		//add any new arcs from the chain of nodes
		if seq == nil || seq.colour != 0 {
			continue
		}
		g.colour(seq) //this will colour a connected component of the graph

		//test for any merging
		for i := 1; i < len(seq.nodes); i++ {
			a := seq.nodes[i-1]
			b := seq.nodes[i]
			if a.offset + a.length > b.offset {
				if a.rc {
					g.mergeNodes(b.node,a.node)
				} else {
					g.mergeNodes(a.node,b.node)

				}
			}
		}
	}

}

func (g *OverlapGraph) colour(seq *SequenceNode) {
	firstColour := g.nextColour
	secondColour := rcColour(g.nextColour)
	if isRCColour(firstColour) {
		secondColour,firstColour = firstColour,secondColour
	}
	if firstColour < secondColour {
		g.nextColour = secondColour+1
	} else {
		g.nextColour = firstColour+1
	}

	open := make([]*SequenceNode, 1, 100)
	open[0] = seq
	seq.colour = firstColour

	tempA := util.NewIntSet()
	tempB := util.NewIntSet()
	tempA2 := util.NewIntSet()
	tempB2 := util.NewIntSet()

	for len(open) > 0 {
		next := open[len(open)-1]
		open = open[:len(open)-1]
		if len(next.nodes) < 2 {
			continue
		}
		otherColour := firstColour
		if next.colour == firstColour {
			otherColour = secondColour
		}
		//find the first coloured node (if any). Link out from there.
		firstNode := -1
		for i, n := range next.nodes {
			if n.node.colour != 0 {
				firstNode = i
				break
			}
		}
		if firstNode == -1 && next != seq {
			continue //this sequence no longer connects coloured nodes
		}
		if firstNode == -1 {
			firstNode = len(next.nodes)-1 //this is just the root sequence case
			fmt.Println("NEW ROOT SEQUENCE")
		}
		fmt.Println("Arcs based on seq ",next.id," chain length ",len(next.nodes))
		before := len(next.nodes)
		//run back from the first node, then up from it
		for i := firstNode; i >= 0; i-- { //all nodes belonging to this sequence
			arc := next.nodes[i]
			newNode := arc.node.colour == 0
			//add arcs up the chain
			if i > 0 && !arc.approximate && !arc.node.isAdjacent(next.nodes[i-1].node) && !next.nodes[i-1].node.isAdjacent(arc.node) {
				g.addArcBetween(next.nodes[i-1].node,arc.node,tempA,tempB,tempA2,tempB2,firstColour)
				//note that the above function can remove arcs from nodes' sequences lists. Including the current node.
			}
			if newNode && arc.node.colour != 0 { //first visit
				for _, childArc := range arc.node.sequences {
					//TODO: don't set child colours? base them on the nodes..
					if !childArc.approximate && childArc.sequence.colour == 0 {
						if arc.rc == childArc.rc { //we flipped back to the parent (next) colour
							childArc.sequence.colour = next.colour
						} else {
							childArc.sequence.colour = otherColour
						}
						open = append(open, childArc.sequence)
					}
				}
			}
		}
		firstNode -= (before-len(next.nodes)) //update, based on deletes
		//more carefully this time: nodes might get removed part way through
		for i := firstNode; i < len(next.nodes); i++ { //all nodes belonging to this sequence
			before := len(next.nodes)
			arc := next.nodes[i]
			newNode := arc.node.colour == 0
			//add arcs up the chain
			if i > 0 && !arc.approximate && !arc.node.isAdjacent(next.nodes[i-1].node) && !next.nodes[i-1].node.isAdjacent(arc.node) {
				g.addArcBetween(next.nodes[i-1].node,arc.node,tempA,tempB,tempA2,tempB2,firstColour)
				//note that the above function can remove arcs from nodes' sequences lists. Including the current node.
			}
			if newNode && arc.node.colour != 0 { //first visit
				for _, childArc := range arc.node.sequences {
					//TODO: don't set child colours? base them on the nodes..
					if !childArc.approximate && childArc.sequence.colour == 0 {
						if arc.rc == childArc.rc { //we flipped back to the parent (next) colour
							childArc.sequence.colour = next.colour
						} else {
							childArc.sequence.colour = otherColour
						}
						open = append(open, childArc.sequence)
					}
				}
			}
			if before != len(next.nodes) { //This node was removed
				i--
			}
		}
		if next == seq && len(open) == 0 {
			//this means our first sequence had no children. Roll back the colours.
			g.nextColour -= 2
		}
	}
}


//GetCoveredSequences tests all sequences and returns whether or not they have overlaps near both ends.
func (g *OverlapGraph) GetCoveredSequences() []bool {
	bc := 0
	fc := 0
	c := 0
	total := 0
	covered := make([]bool, len(g.sequences),len(g.sequences))
	for i, s := range g.sequences {
		if s == nil {
			continue
		}
		if s.covered {
			covered[i] = true
			c++
			bc++
			fc++
		} else if s.coveredFront {
			fc++
		} else if s.coveredBack {
			bc++
		}
		total++
	}
	fmt.Println("Covered:",c,fc,bc,"/",total)
	return covered
}

//Reduces subgraphs to linear chains from the given node to the next convergence
//If no convergence is found by the specified maximum depth, this process halts
/*
func (g *OverlapGraph) linearise(start *Node, maxDepth int, visited []bool) {
	states := NewPriorityQueue()
	states.PushNode(start,0,false)
	currentLength := 0
	layout := make([]int, 0, 20)
	layoutRC := make([]bool, 0, 10)
	for depth := 0; depth < maxDepth && len(states) > 0; depth++ {
		next,distance,rc := states.PopNode()
		fmt.Println("next is",next.id,rc,"at",distance,"with current length",currentLength)
		//check whether the next node already exists
		if visited[next.id] {
			found := false
			for i := len(layout)-1; i > 0; i-= 2 {
				//decide whether or not it is close enough to the old version
				if layout[i] == next.id {
					if rc != layoutRC[i/2] && layout[i-1] + next.Consensus.Len() > distance {
						//it overlaps
						found = true
						break
					}
					fmt.Println("Already visited this one but rc:",rc==layoutRC[i/2],"and old distance:",layout[i-1])
					return //no good. This whole segment is invalid
				}
			}
			if !found {
				fmt.Println("Said visited, but not found in the current layout.")
				//this means we have reached something already linearised.. best stop.
				return
			}
			continue
		}
		visited[next.id] = true
		if distance < currentLength {
			fmt.Println("not far enough past",currentLength)
			return //unmerged overlapping overlaps!
		}
		layout = append(layout,currentLength)
		layout = append(layout,next.id)
		layoutRC = append(layoutRC,rc)
		if len(states) == 0 && len(layout) > 0 {
			//convert current layout's Nodes to a chain
			prevNode := g.nodes[layout[1]]
			prevRC := layoutRC[0]
			length := prevNode.Consensus.Len()
			for i := 3; i < len(layout); i+= 2 {
				n := g.nodes[layout[i]]
				n.in = n.in[:0]
				prevNode.out = prevNode.out[:0]
				g.addArc(prevNode, n, layout[i-1]-length, prevRC, layoutRC[i/2])
				prevRC = layoutRC[i/2]
				length += layout[i-1] + n.Consensus.Len()
			}
			//run forward to the next interesting bit, reset depth
			for len(next.out) == 1 {
				next = next.out[0].to
			}
			layout = layout[:0]
			depth = 0
			currentLength = 0
		}
		currentLength = distance + next.Consensus.Len()
		for _, out := range next.out {
			fmt.Println("Adding child",out.to.id,"at",currentLength,"+",out.length," = ",currentLength+out.length)
			states.PushNode(out.to, currentLength + out.length, rc != (out.fromRC != out.toRC))
		}
	}
}

func (g *OverlapGraph) Linearise() int {
	visited := make([]bool, len(g.nodes), len(g.nodes))
	for _, n := range g.nodes {
		if n != nil && len(n.out) == 1 && len(n.in) != 1 {
			g.linearise(n, 10, visited)
			return n.id
		}
	}
	return -1
}

//GetChains finds contiguous, unambiguous, chains of nodes
func (g *OverlapGraph) GetChains() [][]*Node {
	visited := make([]bool, len(g.nodes), len(g.nodes))
	roots := make([]*Node, 0, 100)
	for i, n := range g.nodes {
		if n != nil {
			if len(n.out) == 1 {
				if len(n.in) != 1 {
					roots = append(roots,n)
				}
			}
			visited[i] = len(n.out) != 1 || len(n.in) != 1
		} else {
			visited[i] = true
		}
	}
	chains := make([][]*Node, 0, len(roots))
	for _, next := range roots {
		//traverse the linear chain
		if visited[next.out[0].to.id] {
			continue //a single node does not a chain make
		}
		chain := make([]*Node, 1, 20)
		chain[0] = next
		next = next.out[0].to
		for !visited[next.id] && len(next.out) != 0 {
			chain = append(chain, next)
			next = next.out[0].to
		}
		chain = append(chain, next)
		chains = append(chains, chain)
	}
	return chains
}
*/
func (g *OverlapGraph) PrintGFA() {
	fmt.Println("H\tVN:Z:1.0")
	for _, n := range g.nodes {
		if n != nil {
			//fmt.Printf("S\t%v\t*\tLN:i:%v\n",n.id,n.Consensus.Len())
			fmt.Printf("S\t%v_%v\t*\tLN:i:%v\n",n.id,n.colour,n.Consensus.Len())
		}
	}
	for _, n := range g.nodes {
		if n == nil {
			continue
		}
		for _, a := range n.out {
			if a.fromRC != a.toRC {
				if a.fromRC {
					fmt.Printf("L\t%v_%v\t-\t%v_%v\t+\t%vM\n",a.from.id,a.from.colour,a.to.id,a.to.colour,a.length)
					//fmt.Printf("L\t%v\t-\t%v\t+\t%vM\n",a.from.id,a.to.id,a.length)
				} else {
					fmt.Printf("L\t%v_%v\t+\t%v_%v\t-\t%vM\n",a.from.id,a.from.colour,a.to.id,a.to.colour,a.length)
					//fmt.Printf("L\t%v\t+\t%v\t-\t%vM\n",a.from.id,a.to.id,a.length)
				}
			} else {
				//fmt.Printf("L\t%v\t+\t%v\t+\t%vM\n",a.from.id,a.to.id,a.length)
				fmt.Printf("L\t%v_%v\t+\t%v_%v\t+\t%vM\n",a.from.id,a.from.colour,a.to.id,a.to.colour,a.length)
			}
		}
	}
}
