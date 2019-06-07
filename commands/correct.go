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

type correctCommand struct {
	args  map[string]string
	alias map[string]string
	desc  map[string]string
}

func NewCorrectCommand() Command {
	args, alias, desc := MakeArgs(
		[]string{"overlap_size", "num_seeds", "seed_batch_size", "chunk_size", "k", "min_hits","num_workers", "input", "trim", "front_adapters", "back_adapters","model","himem"},
		[]string{"1000", "15", "10000", "10000", "10","0.25","4","", "0", "", "","","true"},
		[]string{"Size of overlap to search for in bases", "Minimum number of seeds to generate for each overlap query", "Maximum total unique seeds to use in each query batch", "Size to chop long reads into for querying against, in bases", "Number of bases in each seed", "Minimum proportion of seeds that must match each query", "Number of worker threads to spawn", "Fasta/fastq input file", "Whether to search for and trim adapters: 0=off, 1=on", "Fasta/fastq file containing front adapters", "Fasta/fastq file containing back adapters","K-mer numeric values to use in alignment","Whether to cache all reads in memory"})
	ov := correctCommand{args: args, alias: alias, desc: desc}
	return &ov
}

func (com *correctCommand) GetName() string {
	return "correct"
}

func (com *correctCommand) GetArgs() (map[string]string, map[string]string, map[string]string) {
	return com.args, com.alias, com.desc
}

func (com *correctCommand) Run(args map[string]string) {
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
		numQuerySeqs := 1
		//find the longest reamining sequences
		ids, lengths := seqSet.GetIDsByLength()
		if len(ids) == 0 || lengths[len(lengths)-1] < 1000 {
			break //no more reasonable-length sequences
		}
		//pick some sequences til we think we'll get around the right number of seeds
		last := len(lengths)-1
		start := last
		approxSeeds := (lengths[start]/overlapSize + 1)*numSeeds
		for start >= 0 && approxSeeds < seedBatchSize {
			approxSeeds += (lengths[start]/overlapSize + 1)*numSeeds
			start--
			numQuerySeqs++
		}
		if start < last {
			//TODO: TEMP
			start = last-1
			ids = ids[start+1:]
			fmt.Println("Using",len(lengths)-(start+1),"sequences with lengths",lengths[last],"to",lengths[start+1])
		} else {
			ids = ids[last:]
			fmt.Println("Using one sequence with length",lengths[last])
		}

		fmt.Println("Query ids are ",ids)

		//prepare overlappers, get the next query set
		seqs := seqSet.GetSequencesByID(ids)
		seedIndex := seeds.NewSeedIndex(uint(k))
		overlapper := overlap.NewOverlapper(seedIndex, chunkSize, numWorkers, overlapSize, 10, hitFraction)
		queries := overlapper.PrepareQueries(numSeeds, seedBatchSize, values, seqs, overlap.QueryAll)

		fmt.Println("Produced a query set of", len(queries), "queries using", seedIndex.Size(), "seeds.")

		// 1. get query results (collated by long query, ordered)
		results := performQueries(queries, overlapper, overlapSize, seedIndex, seqSet, ids)
		c1 := 0 //count of all shorter reads found matching a query
		sentCount := 0
		successfulCount := 0
		seedConsensus := make([][]*overlap.SeedContig, len(results))
		seqIds := util.NewIntSetCapacity(seqSet.Size())
		for j, rs := range results {
			for _,hits := range rs {
				//NOTE: this also sorts the hits in a non-useful order!
				removeDuplicates(&hits)
				c1 += len(hits)
			}
			// 2. remove any that don't match across whole length
			sort.Sort(byOffset(rs))
			overlap.CleanupOverlaps(rs, overlapSize, k)

			seedConsensus[j] = seedSpaceConsensus(rs,seedIndex,seqIds)
		}
		//Note: at this point the contigs are seed sequences aligned to a "seed consensus" with id of the original long sequence

		fmt.Println("Attempted",sentCount,"seed-space consensuses, succeeding with",successfulCount)

		//Throw out the seed index. Its job is done.
		seedIndex.Destroy()
		allSeq := getAllSequences(seqIds, seqSet)

		fmt.Println("Preparing base-space consensus of all query results.")
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
			queriesReady <- true
		}()
		//wait for all consensus sequences to be completed before closing the output channel
		for i := 0; i < numWorkers; i++ {
			<-done
		}
		close(allConsensus)

		<-queriesReady
		refinedResults := performQueries(nextQueries, overlapper, overlapSize, seedIndex, seqSet, ids)
		//prepare for seed-space consensus of the second round. Overwrite with new allocations to help garbage collection
		seedConsensus = make([][]*overlap.SeedContig, len(refinedResults))
		seqIds.Clear()
		fmt.Println("FINAL query results (",len(refinedResults),"):")
		for j, rs := range refinedResults {
			fmt.Println("Refined results ",j," has ",len(rs),"contigs")
			for _,hits:= range rs {
				removeDuplicates(&hits)
			}
			fmt.Println("After removing duplicates:",len(rs),"contigs")
			sort.Sort(byOffset(rs)) //not sure whether we need this here. Can't hurt.
			overlap.CleanupOverlaps(rs, overlapSize, k)
			seedConsensus[j] = seedSpaceConsensus(rs,seedIndex,seqIds)
			/*for i, hits := range rs {
				fmt.Print(i)
				for _,m := range hits {
					fmt.Print(",",m.SeqB.GetID())
				}
				fmt.Println()
			}*/
		}
		seedIndex.Destroy()
		// 4. Generate pileup of "non-flapping" sequences, mark flapping regions
		overlap.NewPileup(seedConsensus[0])

		// 5. Take consensus over longer contiguous regions
		// 6. Flag to ignore short reads within non-flapping regions
		// 7. Output the corrected reads

		//"pileup" package, or go in overlap/pileup.go <-
		// handle the "series of subsequences" format
		// - clean up those that aren't in order (from graph.go)
		// - test for overhangs (from graph.go)
		// - select subsequences for max length contiguous segments with coverage @~90% max

		//TODO: move the covered test out of the overlap graph
		// as a per-sequence check
		/*covered := overlapGraph.GetCoveredSequences()
		coveredCount := 0
		for i := 0; i < len(covered); i++ {
			if covered[i] {
				seqSet.SetIgnore(i,true)
				coveredCount++
			}
		}*/
		//fmt.Println("Ignoring ", coveredCount, " of the future sequences.")
		break
	}

	//fmt.Println("tried to linearise from:",overlapGraph.Linearise())
	pprof.StopCPUProfile()
	f.Close()
}

func seedSpaceConsensus(rs [][]*seeds.SeedMatch, seedIndex *seeds.SeedIndex, seqIds *util.IntSet) []*overlap.SeedContig {
	seedConsensus := make([]*overlap.SeedContig, len(rs))
	for i, hits := range rs {
		if len(hits) >= 3 {
			contig := overlap.BuildConsensus(seedIndex, hits)
			if contig != nil && len(contig.Parts) >= 3 {
				seedConsensus[i] = contig
				//take note of the required sequences for later base-space consensus
				for _, part := range contig.Parts {
					seqIds.Add(uint(part))
				}
				//slight hack: give the seed consensus the id of the original query
				originalID := hits[0].SeqA.GetID()
				contig.Combined.SetID(originalID)
				//which hit is on the original?
				original := -1
				for k, part := range contig.Parts {
					if part == originalID {
						original = k
						break
					}
				}
				//and we also need to set its offsets
				if original == -1 {
					fmt.Println("WARNING: seed consensus does not contain original sequence.")
					contig.Combined.SetOffsets(hits[0].SeqA.GetOffset(), hits[0].SeqA.GetInset())
				} else {
					//TODO: recalculate the inset too
					contig.Combined.SetOffsets(hits[0].SeqA.GetOffset() + contig.Offsets[original], hits[0].SeqA.GetInset())
				}
			}
		}
	}
	return seedConsensus
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
	//run through the queries, they will be collated by sequence, we'll put them in order later
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
//sort a list of list of matches by their shared target sequence
type byOffset [][]*seeds.SeedMatch

func (a byOffset) Len() int           { return len(a) }
func (a byOffset) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a byOffset) Less(i, j int) bool {
	if len(a[i]) == 0 {
		return len(a[j]) == 0
	} else if len(a[j]) == 0 {
		return false
	}
	return a[i][0].SeqA.GetOffset() < a[j][0].SeqA.GetOffset()
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
				if i < n {
					copy(rs[i:],rs[i+1:])
				}
				//rs[i] = rs[n]
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
		if _, cons := consensus.BuildConsensus(contig, sequences, model, false); cons != nil {
			//fmt.Println(cons.String())
			if output != nil {
				output <- cons
			}
		}
	}
	done <- true
}
