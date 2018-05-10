package commands

import (
	"fmt"
	"github.com/jteutenberg/downpore/consensus"
	"github.com/jteutenberg/downpore/model"
	"github.com/jteutenberg/downpore/overlap"
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/trim"
	"os"
	"runtime/pprof"
	"sort"
)

type fullmonteCommand struct {
	args  map[string]string
	alias map[string]string
	desc  map[string]string
}

func NewFullMonteCommand() Command {
	args, alias, desc := MakeArgs(
		[]string{"overlap_size", "num_seeds", "seed_batch_size", "chunk_size", "k", "min_hits","num_workers", "input", "trim", "front_adapters", "back_adapters","model","himem"},
		[]string{"1000", "15", "10000", "10000", "10","0.25","4","", "0", "", "","data/model.txt","true"},
		[]string{"Size of overlap to search for in bases", "Minimum number of seeds to generate for each overlap query", "Maximum total unique seeds to use in each query batch", "Size to chop long reads into for querying against, in bases", "Number of bases in each seed", "Minimum proportion of seeds that must match each query", "Number of worker threads to spawn", "Fasta/fastq input file", "Whether to search for and trim adapters: 0=off, 1=on", "Fasta/fastq file containing front adapters", "Fasta/fastq file containing back adapters","K-mer numeric values to use in alignment","Whether to cache all reads in memory"})
	ov := fullmonteCommand{args: args, alias: alias, desc: desc}
	return &ov
}

func (com *fullmonteCommand) GetName() string {
	return "fullmonte"
}

func (com *fullmonteCommand) GetArgs() (map[string]string, map[string]string, map[string]string) {
	return com.args, com.alias, com.desc
}

func (com *fullmonteCommand) Run(args map[string]string) {
	overlapSize := ParseInt(args["overlap_size"])
	numSeeds := ParseInt(args["num_seeds"]) //minimum seeds per overlap region
	seedBatchSize := int(ParseInt(args["seed_batch_size"]))
	chunkSize := uint(ParseInt(args["chunk_size"])) //slice reads into chunks of this size
	k := ParseInt(args["k"])
	hitFraction := ParseFloat(args["min_hits"])
	numWorkers := ParseInt(args["num_workers"])
	model := model.NewModel(args["model"], false)

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
		end := 1
		approxSeeds += (lengths[0]/overlapSize + 1)*numSeeds
		for end < len(ids) && approxSeeds < seedBatchSize {
			approxSeeds += (lengths[end]/overlapSize + 1)*numSeeds
			end++
			numQuerySeqs++
		}
		if end > 1 {
			ids = ids[:end-1]
			fmt.Println("Using",end,"sequences with lengths",lengths[0],"to",lengths[end-1])
		} else {
			ids = ids[:1]
			fmt.Println("Using one sequence with length",lengths[0])
		}

		//prepare overlappers, get the next query set
		seqs := seqSet.GetSequencesByID(ids)
		seedIndex := seeds.NewSeedIndex(uint(k))
		overlapper := overlap.NewOverlapper(seedIndex, chunkSize, numWorkers, overlapSize, 10, hitFraction)
		queries := overlapper.PrepareQueries(numSeeds, seedBatchSize, values, seqs, overlap.QueryAll) //TODO: separate out query construction so a custom one can be used

		fmt.Println("Produced a query set of", len(queries), "queries using", seedIndex.Size(), "seeds.")

		//TODO:
		// - for each (long?) query sequence, lay out the hits over time
		// - look for any that flap. Mark the attachment as potential repeat location. Discard the flaps? Keep, but put in separate contig, index the repeat?
		// - on remaining good sequences, perform consensus
		// - re-query?
		// - remove contained sequences from future queries

		//now generate the index of all sequences using these seeds
		contigs := overlapsFromQueries(queries, overlapper, seedIndex, seqSet, numQuerySeqs)
		for i,contig := range contigs {
			fmt.Println(i,contig.Parts)
		}
		//genomeSize = seqSet.GetBases() / int64(coverage)

		if overlapGraph == nil {
			//TODO: makes more sense to get the number of sequences from the sequence set
			overlapGraph = overlap.NewOverlapGraph(int(seedIndex.GetNumSequences()))
		}

		//Throw out the seed index. Its job is done.
		seedIndex.Destroy()

		allSeq := getAllSequences(contigs, seqSet)

		fmt.Println("Preparing consensus of all query results.")
		allConsensus := make(chan sequence.Sequence, numWorkers)
		contigsIn := make(chan *overlap.SeedContig, numWorkers)
		//setup workers to prepare consensus sequences from overlaps
		done := make(chan bool, numWorkers+1)
		for i := 0; i < numWorkers; i++ {
			go consensusWorker(contigsIn, allSeq, model, done, allConsensus, nil)
		}
		//send overlaps in to have their consensus made
		go func() {
			for _, contig := range contigs {
				contigsIn <- contig
			}
			close(contigsIn)
		}()
		//wait for all consensus sequences to be completed before closing the output channel
		go func() {
			for i := 0; i < numWorkers; i++ {
				<-done
			}
			close(allConsensus)
		}()

		//prepare new queries from the consensus sequences
		/*seedIndex = seeds.NewSeedIndex(uint(k))
		overlapper = overlap.NewOverlapper(seedIndex, chunkSize, numWorkers, overlapSize*10, 10, hitFraction) //oversized overlaps: don't split queries
		queries = overlapper.PrepareQueries(numSeeds, math.MaxInt32, values, allConsensus, false) //no limit on number of seeds

		fmt.Println("Round 2 will query using", len(queries), "from", len(contigs), "potential consensus sequences.")

		//at this point the queries have drained allConsensus and are ready to be used

		contigs = overlapsFromQueries(queries, overlapper, seedIndex, seqSet, numQuerySeqs, coverage)

		//take the consensus again (hopefully improved) and add to the graph
		allSeq = getAllSequences(contigs, seqSet)
		contigsIn = make(chan *overlap.SeedContig, numWorkers)
		//setup workers to prepare consensus sequences from overlaps
		done = make(chan bool, numWorkers+1)
		for i := 0; i < numWorkers; i++ {
			go consensusWorker(contigsIn, allSeq, model, done, nil, overlapGraph)
		}
		//send overlaps in to have their consensus made
		go func() {
			for _, contig := range contigs {
				contigsIn <- contig
			}
			close(contigsIn)
		}()
		//wait for everything to be added to the graph
		for i := 0; i < numWorkers; i++ {
			<-done
		}

		bridges := overlapGraph.GetBridgableContigs(uint(20))
		fmt.Println("Found ", len(bridges), "bridgable locations")
		/*contigsIn = make(chan *overlap.SeedContig, numWorkers*2)
		done = make(chan bool, numWorkers+1)
		for i := 0; i < numWorkers; i++ {
			go bridgeWorker(contigsIn, allSeq, model, overlapGraph, done)
		}
		go func() {
			for _, contig := range bridges {
				contigsIn <- contig
			}
			close(contigsIn)
		}()
		for i := 0; i < numWorkers; i++ {
			<-done
		}*/

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

func overlapsFromQueries(queries []*overlap.SeedQuery, overlapper overlap.Overlapper, seedIndex *seeds.SeedIndex, seqSet sequence.SequenceSet, numQuerySeqs int) []*overlap.SeedContig {
	seqs := seqSet.GetSequences()
	overlapper.AddSequences(seqs)
	fmt.Println("Finished indexing.")

	//Do all the queries, clear the index
	queryResults := make([][]*seeds.SeedMatch, 1, numQuerySeqs)
	matches := overlapper.FindOverlaps(queries)
	for match := range matches {
		qId := match.QueryID
		for qId >= len(queryResults) {
			queryResults = append(queryResults, nil)
		}
		if queryResults[qId] == nil {
			queryResults[qId] = make([]*seeds.SeedMatch, 0, 20)
		}
		//this collates forward and reverse-complement queries
		queryResults[qId] = append(queryResults[qId], match)
	}
	fmt.Println("Query completed. ",len(queryResults),"total queries used.")
	//each query has now found all matching sequences. Form seed-consensus.
	contigs := make([]*overlap.SeedContig, 0, len(queryResults))
	for _, results := range queryResults {
		if results != nil && len(results) > 1 {
			//remove any duplicate results (happens when two overlapping read chunks get hit, or local repeats)
			removeDuplicates(&results)
			//seed-space consensus
			contig := overlap.BuildConsensus(seedIndex, results)
			if contig != nil && len(contig.Parts) > 2 {
				contigs = append(contigs, contig)
			}
		}
	}
	queryResults = queryResults[:0]
	return contigs
}

//a sortable version of seed match slices
type byID []*seeds.SeedMatch

func (a byID) Len() int           { return len(a) }
func (a byID) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a byID) Less(i, j int) bool { return a[i].SeqB.GetID() < a[j].SeqB.GetID() }

func removeDuplicates(results *[]*seeds.SeedMatch) {
	rs := *results
	//sort the slice
	sort.Sort(byID(rs))
	//check for duplicates and remove
	prev := rs[len(rs)-1]
	for i := len(rs) - 2; i >= 0; i-- {
		m := rs[i]
		if m.SeqB.GetID() == prev.SeqB.GetID() {
			n := len(rs) - 1
			rs[i] = rs[n]
			rs = rs[:n]
		}
		prev = m
	}
	*results = rs
}

func getAllSequences(contigs []*overlap.SeedContig, sequences sequence.SequenceSet) []sequence.Sequence {
	ids := make([]bool, sequences.Size(), sequences.Size())
	count := 0
	for _, contig := range contigs {
		for _, id := range contig.Parts {
			if !ids[id] {
				ids[id] = true
				count++
			}
		}
	}
	idList := make([]int, 0, len(ids))
	for id, exists := range ids {
		if exists {
			idList = append(idList, id)
		}
	}
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

func consensusWorker(contigs <-chan *overlap.SeedContig, sequences []sequence.Sequence, model model.Model, done chan<- bool, output chan<- sequence.Sequence, graph *overlap.OverlapGraph) {
	for contig := range contigs {
		if _, cons := consensus.BuildConsensus(contig, sequences, model, false); cons != nil {
			if output != nil {
				output <- cons
			}
			if graph != nil {
				graph.AddNode(contig, cons)
			}
		}
	}
	done <- true
}
func bridgeWorker(contigs <-chan *overlap.SeedContig, sequences []sequence.Sequence, model model.Model, graph *overlap.OverlapGraph, done chan<- bool) {
	for contig := range contigs {
		if contig, cons := consensus.BuildConsensus(contig, sequences, model, true); cons != nil {
			graph.AddBridge(contig, cons)
		}
	}
	done <- true
}
