package commands

import (
	"fmt"
	"github.com/jteutenberg/downpore/consensus"
	"github.com/jteutenberg/downpore/model"
	"github.com/jteutenberg/downpore/overlap"
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/trim"
	"github.com/jteutenberg/downpore/util"
	"os"
	"runtime/pprof"
	"sort"
)

type fulldenovoCommand struct {
	args  map[string]string
	alias map[string]string
	desc  map[string]string
}

func NewFullDenovoCommand() Command {
	args, alias, desc := MakeArgs(
		[]string{"overlap_size", "num_seeds", "seed_batch_size", "chunk_size", "k", "min_hits","num_workers", "input", "trim", "front_adapters", "back_adapters","model","himem"},
		[]string{"1000", "15", "10000", "10000", "10","0.25","4","", "0", "", "","","true"},
		[]string{"Size of overlap to search for in bases", "Minimum number of seeds to generate for each overlap query", "Maximum total unique seeds to use in each query batch", "Size to chop long reads into for querying against, in bases", "Number of bases in each seed", "Minimum proportion of seeds that must match each query", "Number of worker threads to spawn", "Fasta/fastq input file", "Whether to search for and trim adapters: 0=off, 1=on", "Fasta/fastq file containing front adapters", "Fasta/fastq file containing back adapters","K-mer numeric values to use in alignment","Whether to cache all reads in memory"})
	ov := fulldenovoCommand{args: args, alias: alias, desc: desc}
	return &ov
}

func (com *fulldenovoCommand) GetName() string {
	return "fulldenovo"
}

func (com *fulldenovoCommand) GetArgs() (map[string]string, map[string]string, map[string]string) {
	return com.args, com.alias, com.desc
}

func (com *fulldenovoCommand) Run(args map[string]string) {
	overlapSize := ParseInt(args["overlap_size"])
	numSeeds := ParseInt(args["num_seeds"]) //minimum seeds per overlap region
	seedBatchSize := int(ParseInt(args["seed_batch_size"]))
	chunkSize := uint(ParseInt(args["chunk_size"])) //slice reads into chunks of this size
	k := ParseInt(args["k"])
	hitFraction := ParseFloat(args["min_hits"])
	numWorkers := ParseInt(args["num_workers"])
	var mod model.Model
	if args["model"] != "" {
		mod = model.NewModel(args["model"], false)
	}

	var overlapGraph *overlap.OverlapGraph

	seqSet := sequence.NewFastaSequenceSet(args["input"], overlapSize, numWorkers, ParseBool(args["himem"]), false)

	if args["trim"] == "1" {
		trimmer := trim.LoadTrimmer(args["front_adapters"], args["back_adapters"], 5)
		trimmer.Trim(seqSet, numWorkers)
		trimmer.PrintStats(seqSet)
	}

	values := getKmerValues(args["seed_values"],k,numWorkers, seqSet)

	f, _ := os.Create("./cprof")
	pprof.StartCPUProfile(f)
	for ;; {
		//first select any extensions from the graph?
		approxSeeds := 0
		numQuerySeqs := 1
		//then select some long sequences to top up with
		ids, lengths := seqSet.GetIDsByLength()
		if len(ids) == 0 {
			break //no more sequences
		}
		last := len(lengths)-1
		start := last
		approxSeeds += (lengths[start]/overlapSize + 1)*numSeeds
		for start >= 0 && approxSeeds < seedBatchSize {
			approxSeeds += (lengths[start]/overlapSize + 1)*numSeeds
			start--
			numQuerySeqs++
		}
		if start < last {
			ids = ids[start+1:]
			fmt.Println("Using",len(lengths)-start,"sequences with lengths",lengths[last],"to",lengths[start+1])
		} else {
			ids = ids[last:]
			fmt.Println("Using one sequence with length",lengths[last])
		}

		//prepare overlappers, get the next query set
		seqs := seqSet.GetSequencesByID(ids)
		seedIndex := seeds.NewSeedIndex(uint(k))
		overlapper := overlap.NewOverlapper(seedIndex, chunkSize, numWorkers, overlapSize, 10, hitFraction)
		queries := overlapper.PrepareQueries(numSeeds, seedBatchSize, values, seqs, overlap.QueryAll)

		fmt.Println("Produced a query set of", len(queries), "queries using", seedIndex.Size(), "seeds.")


		// 1. get query results (collated by long query, ordered)
		results := performQueries(queries, overlapper, overlapSize, seedIndex, seqSet, ids)
		c1 := 0
		c2 := 0
		sentCount := 0
		successfulCount := 0
		seedConsensus := make([][]*overlap.SeedContig, len(results))
		seqIds := util.NewIntSetCapacity(seqSet.Size())
		for j, rs := range results {
			for _,hits:= range rs {
				removeDuplicates(&hits)
				c1 += len(hits)
			}
			// 2. remove any that don't match across whole length
			overlap.CleanupOverlaps(rs, overlapSize, k)
			seedConsensus[j] = make([]*overlap.SeedContig, len(rs))
			for i, hits := range rs {
				fmt.Print(i)
				for _,m := range hits {
					fmt.Print(",",m.SeqB.GetID())
				}
				fmt.Println()
				c2 += len(hits)
				// 3. consensus, left to right
				// a. seed-space consensus first
				if len(hits) >= 3 {
					sentCount++
					contig := overlap.BuildConsensus(seedIndex, hits)
					if contig != nil && len(contig.Parts) >= 3 {
						successfulCount++
						seedConsensus[j][i] = contig
						//take note of the required sequences for later base-space consensus
						for _, part := range contig.Parts {
							seqIds.Add(uint(part))
						}
					}
				}
			}
		}
		fmt.Println("Attempted ",sentCount,"consensuses, succeeding with",successfulCount)

		if overlapGraph == nil {
			//TODO: makes more sense to get the number of sequences from the sequence set
			overlapGraph = overlap.NewOverlapGraph(int(seedIndex.GetNumSequences()))
		}

		//Throw out the seed index. Its job is done.
		seedIndex.Destroy()
		allSeq := getAllSequences(seqIds, seqSet)

		fmt.Println("Preparing consensus of all query results.")
		allConsensus := make(chan sequence.Sequence, numWorkers)
		contigsIn := make(chan *overlap.SeedContig, numWorkers)
		//setup workers to prepare consensus sequences from overlaps
		done := make(chan bool, numWorkers+1)
		for i := 0; i < numWorkers; i++ {
			go consensusWorker(contigsIn, allSeq, mod, done, allConsensus)
		}
		//send overlaps in to have their consensus made
		go func() {
			sentCount := 0
			for _, contigs := range seedConsensus {
				for _, contig := range contigs {
					if contig != nil {
						sentCount++
						contigsIn <- contig
					}
				}
			}
			fmt.Println("Sent through",sentCount,"overlaps for consensus building")
			close(contigsIn)
		}()

		// Get seeds from the consensus sequences to use as a next round of queries
		seedIndex = seeds.NewSeedIndex(uint(k))
		overlapper = overlap.NewOverlapper(seedIndex, chunkSize, numWorkers, overlapSize, 10, hitFraction)

		queryIn := make(chan sequence.Sequence, numWorkers*2)
		qCount := 0
		go func() {
			for c := range allConsensus {
				queryIn <- c
				qCount++
			}
			fmt.Println("Received",qCount,"consensus results.")
			close(queryIn)
		}()
		//and read the queries out
		queriesReady := make(chan bool, 1)
		var nextQueries []*overlap.SeedQuery
		go func() {
			nextQueries = overlapper.PrepareQueries(numSeeds, seedBatchSize, values, queryIn, overlap.QueryAll)
			fmt.Println("Turned",len(nextQueries)/2,"consensus sequences into queries.")
			queriesReady <- true
		}()
		//wait for all consensus sequences to be completed before closing the output channel
		for i := 0; i < numWorkers; i++ {
			<-done
		}
		close(allConsensus)

		<-queriesReady
		refinedResults := performQueries(nextQueries, overlapper, overlapSize, seedIndex, seqSet, ids)
		seedIndex.Destroy()
		fmt.Println("FINAL:")
		for _, rs := range refinedResults {
			for _,hits:= range rs {
				removeDuplicates(&hits)
			}
			overlap.CleanupOverlaps(rs, overlapSize, k)
			for i, hits := range rs {
				fmt.Print(i)
				for _,m := range hits {
					fmt.Print(",",m.SeqB.GetID())
				}
				fmt.Println()
			}
		}
		//TODO:
		// . fill in any missing sequences to get full sets
		// 5. re-query
		// 6. repeat checks but save flaps etc. Goes into the graph

		//prepare new queries from the consensus sequences

		//turn off any future sequences that are in the graph
		covered := overlapGraph.GetCoveredSequences()
		coveredCount := 0
		for i := 0; i < len(covered); i++ {
			if covered[i] {
				seqSet.SetIgnore(i,true)
				coveredCount++
			}
		}
		fmt.Println("Ignoring ", coveredCount, " of the future sequences.")
		break
	}
	overlapGraph.GenerateArcs()
	overlapGraph.PrintGFA()

	//fmt.Println("tried to linearise from:",overlapGraph.Linearise())
	pprof.StopCPUProfile()
	f.Close()
}

//performQueries indexes all sequences and finds matches for all queries. Matches are returned in in-order collections per query sequence.
//The results are slices of form [query sequence][overlap][hits]
func performQueries(queries []*overlap.SeedQuery, overlapper overlap.Overlapper, overlapSize int, seedIndex *seeds.SeedIndex, seqSet sequence.SequenceSet, querySequences []int) [][][]*seeds.SeedMatch {
	seqs := seqSet.GetSequences()
	overlapper.AddSequences(seqs)

	queryResults := make([][][]*seeds.SeedMatch, len(querySequences))
	queryIndices := make(map[int]int)
	for i, id := range querySequences {
		queryResults[i] = make([][]*seeds.SeedMatch, 0, seqSet.GetLength(id)/overlapSize + 1)
	}
	//because the queries are in-order, we can run through them..
	index := 0
	prevSeq := -1
	for _, q := range queries {
		if q.SequenceID != prevSeq {
			prevSeq = q.SequenceID
			index = 0
		}
		queryIndices[q.ID] = index/2 //query pairs (forward + reverse-complement)
		index++
	}
	//find the overlap matches and assign them to their sequence+index
	matches := overlapper.FindOverlaps(queries)
	for match := range matches {
		seqId := match.SeqA.GetID() //the query sequence
		seqIndex := 0
		for i,index := range querySequences {
			if seqId == index {
				seqIndex = i
				break
			}
		}
		index = queryIndices[match.QueryID]
		for len(queryResults[seqIndex]) <= index {
			queryResults[seqIndex] = append(queryResults[seqIndex], make([]*seeds.SeedMatch,0,20))
		}
		//fmt.Println("Result to seq ",seqIndex,"/",len(queryResults)," at",index,"/",len(queryResults[seqIndex]))
		queryResults[seqIndex][index] = append(queryResults[seqIndex][index], match)
	}
	return queryResults
}

/*
//refineQueries queries each consensus sequence is queried and the results replace entries in oldResults. 
func refineQueries(queries []*overlap.SeedQuery, overlapper overlap.Overlapper, overlapSize int, seedIndex *seeds.SeedIndex, seqSet sequence.SequenceSet, querySequences []int) [][][]*seeds.SeedMatch {
func refineQueries(consensus [][]sequence.Sequence, oldResults [][][]*seeds.SeedMatch, overlapSize int, seqSet sequence.SequenceSet, numWorkers int) {

	//2. Index ...
	seqs := seqSet.GetSequences()
	overlapper.AddSequences(seqs)

	//3. ... and perform the query
	matches := overlapper.FindOverlaps(queries)
	for match := range matches {
		seqId := match.SeqA.GetID() //the query sequence
		fmt.Println("Query ",seqId)
		for i,p := range match.Parts {
			fmt.Println(i,":",p)
		}
	}
}*/

//a sortable version of seed match slices
type byID []*seeds.SeedMatch

func (a byID) Len() int           { return len(a) }
func (a byID) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a byID) Less(i, j int) bool {
	idA := a[i].SeqB.GetID()
	idB := a[j].SeqB.GetID()
	if idA == idB {
		//sort by offset
		return a[i].SeqB.GetOffset() < a[j].SeqB.GetOffset()
	}
	return idA < idB
}

func removeDuplicates(results *[]*seeds.SeedMatch) {
	rs := *results
	//sort the slice
	sort.Sort(byID(rs))
	//check for duplicates and remove
	prev := rs[len(rs)-1]
	for i := len(rs) - 2; i >= 0; i-- {
		m := rs[i]
		if m.SeqB.GetID() == prev.SeqB.GetID() {
			//need an additional check: roughly same part of the query, SeqB?
			centre1 := (m.SeqB.GetOffset() + m.SeqB.GetLength())/2
			centre2 := (prev.SeqB.GetOffset() + prev.SeqB.GetLength())/2
			if (centre1 > prev.SeqB.GetOffset() && centre1-prev.SeqB.GetOffset() < prev.SeqB.GetLength()) || (centre2 > m.SeqB.GetOffset() && centre2-m.SeqB.GetOffset() < m.SeqB.GetLength()) {
				n := len(rs) - 1
				rs[i] = rs[n]
				rs = rs[:n]
			}
		}
		prev = m
	}
	*results = rs
}

func getAllSequences(ids *util.IntSet, sequences sequence.SequenceSet) []sequence.Sequence {
	idList := ids.AsInts()
	if len(idList) == 0 {
		return make([]sequence.Sequence, 0, 0)
	}
	allSeq := make([]sequence.Sequence, idList[len(idList)-1]+1, idList[len(idList)-1]+1)
	seqs := sequences.GetSequencesByID(idList)
	for s := range seqs {
		allSeq[s.GetID()] = s
	}
	return allSeq
}

func consensusWorker(contigs <-chan *overlap.SeedContig, sequences []sequence.Sequence, model model.Model, done chan<- bool, output chan<- sequence.Sequence) {
	for contig := range contigs {
		/*for _, m := range contig.Matches {
			fmt.Println(m.MatchA, m.MatchB, m.SeqB.GetSeedOffset(m.MatchB[0],10)) // consensus seeds matched
			fmt.Println(m.SeqB.String())
		}
		for i, part := range contig.Parts {
			end := contig.Offsets[i]+contig.Lengths[i]
			start := contig.Offsets[i]
			if end > sequences[part].Len() {
				fmt.Println("over length")
				if end > sequences[part].Len()+10 {
					continue
				}
				end = sequences[part].Len()
			}
			if start < 0 {
				fmt.Println("under start")
				if start < -10 {
					continue
				}
				start = 0
			}
			s := sequences[part].SubSequence(start,end)
			if contig.ReverseComplement[i] {
				s = s.ReverseComplement()
			}
			fmt.Println("@",i,contig.ReverseComplement[i])
				fmt.Println(s.String())
			fmt.Println("+")
			qs := s.Quality()
			if qs != nil {
				for j,q := range qs {
					qs[j] = q+33
				}
				fmt.Println(string(qs))
			} else {
				for j := 0; j < s.Len(); j++ {
					fmt.Print('A')
				}
				fmt.Println()
			}
		}*/
		if _, cons := consensus.BuildConsensus(contig, sequences, model, false); cons != nil {
			//fmt.Println(cons.String())
			if output != nil {
				output <- cons
				//TODO: we want the start/end positions for each component sequence too
			}
		}
	}
	done <- true
}
