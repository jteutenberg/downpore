package commands

import (
	"fmt"
	"github.com/jteutenberg/downpore/overlap"
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/util/sequtil"
	"log"
	"math"
	"sync"
	"time"
	"os"
	"runtime/pprof"
)

type overlapCommand struct {
	args  map[string]string
	alias map[string]string
	desc  map[string]string
}

func NewOverlapCommand() Command {
	args, alias, desc := MakeArgs(
		[]string{"overlap_size", "k", "num_seeds", "seed_batch_size", "chunk_size", "num_workers", "input", "himem"},
		[]string{"1000", "10", "15", "10000", "10000", "4", "", "true"},
		[]string{"Size of overlap to search for in bases", "Number of bases in each seed", "Minimum number of seeds to generate for each overlap query", "Maximum total unique seeds to use in each query batch", "Size to chop long reads into for querying against, in bases", "Number of worker threads to spawn", "Fasta/fastq input file", "Whether to cache all reads in memory"})
	ov := overlapCommand{args: args, alias: alias, desc: desc}
	return &ov
}

func (com *overlapCommand) GetName() string {
	return "overlap"
}

func (com *overlapCommand) GetArgs() (map[string]string, map[string]string, map[string]string) {
	return com.args, com.alias, com.desc
}

func (com *overlapCommand) Run(args map[string]string) {
	overlapSize := ParseInt(args["overlap_size"])
	numSeeds := ParseInt(args["num_seeds"]) //minimum seeds per overlap region
	seedBatchSize := int(ParseInt(args["seed_batch_size"]))
	chunkSize := uint(ParseInt(args["chunk_size"])) //slice reads into chunks of this size
	numWorkers := ParseInt(args["num_workers"])
	k := ParseInt(args["k"])

	//run through the sequences file, generating batches of queries

	seqSet := sequence.NewFastaSequenceSet(args["input"], overlapSize, numWorkers, ParseBool(args["himem"]), true)

	log.Printf("Counting all %v-mers in the input...\n", k)
	kmerCounts := sequtil.KmerOccurrences(seqSet.GetSequences(), k, numWorkers)
	values := make([]float64, len(kmerCounts))
	for i, count := range kmerCounts {
		j := seeds.ReverseComplement(uint(i), uint(k))
		ratio := float64(count) / float64(kmerCounts[j]+1)
		if ratio < 1.0 {
			ratio = 1.0 / ratio
		}
		values[i] = ratio //replacing values
	}

	bottom, top := sequtil.TopOccurrences(kmerCounts, uint(k), len(kmerCounts)/100, len(kmerCounts)/50)
	for _, x := range bottom {
		values[x] = 10000
	}
	for _, x := range top {
		values[x] = 10000
	}
	log.Println("Counting complete. Starting indexing and querying...")

	f, _ := os.Create("./cprof")
	pprof.StartCPUProfile(f)

	firstSequence := 0
	start := time.Now()
	for round := 0; ; round++ {
		if round == 1 {
			roundTime := time.Since(start)
			numRounds := float64(seqSet.Size()) / float64(firstSequence)
			if numRounds > 1.5 {
				log.Println("Expected remaining execution time:", roundTime*(time.Duration(int(numRounds+0.5))))
			}
		}
		var seedIndex *seeds.SeedIndex
		seqs := seqSet.GetNSequencesFrom(firstSequence, int(math.MaxInt32))
		var overlapper overlap.Overlapper
		seedIndex = seeds.NewSeedIndex(uint(k))
		overlapper = overlap.NewOverlapper(seedIndex, chunkSize, numWorkers, overlapSize, numSeeds, 0.25)

		queries := overlapper.PrepareQueries(numSeeds, seedBatchSize, values, seqs, false)
		if len(queries) == 0 {
			break
		}

		numQuerySeqs := 0
		firstSequence = queries[len(queries)-1].SequenceID + 1
		for _, q := range queries {
			if q.ID >= numQuerySeqs {
				numQuerySeqs = q.ID + 1
			}
			if q.SequenceID >= firstSequence {
				firstSequence = q.SequenceID + 1
			}
		}

		//now continue through the file, building up a set to query against
		seqs = seqSet.GetSequences()
		queryResults := make([][]*seeds.SeedMatch, numQuerySeqs, numQuerySeqs)

		overlapper.AddSequences(seqs)
		if round == 0 {
			log.Println("Using query sets of around", firstSequence, "sequences against", seqSet.Size(), "sequences.")
		}
		//Do all the queries, clear the index
		matches := overlapper.FindOverlaps(queries)
		hits := 0
		qHits := 0
		for match := range matches {
			hits++
			qId := match.QueryID
			if queryResults[qId] == nil {
				queryResults[qId] = make([]*seeds.SeedMatch, 0, 20)
			}
			if len(queryResults[qId]) == 1 {
				qHits++
			}
			//this collates forward and reverse-complement queries
			queryResults[qId] = append(queryResults[qId], match)
		}

		allResults := make(chan []*seeds.SeedMatch, numWorkers*2)
		done := make(chan bool, numWorkers)
		var lock sync.Mutex
		for i := 0; i < numWorkers; i++ {
			go finalCheckWorker(allResults, seedIndex, seqSet, overlapSize, &lock, done)
		}

		for _, results := range queryResults {
			if results != nil && len(results) > 1 {
				allResults <- results
			}
		}
		close(allResults)
		for i := 0; i < numWorkers; i++ {
			<-done
		}
		queryResults = queryResults[:0]

		seedIndex.Destroy()
	}
	pprof.StopCPUProfile()
	f.Close()
}

func finalCheckWorker(overlaps <-chan []*seeds.SeedMatch, seedIndex *seeds.SeedIndex, seqSet sequence.SequenceSet, overlapSize int, printLock *sync.Mutex, done chan<- bool) {
	k := int(seedIndex.GetSeedLength())
	for results := range overlaps {
		//ensure the seed sequences map to one another as well
		contig := overlap.BuildConsensus(seedIndex, results)
		if contig != nil && len(contig.Parts) > 1 {
			if contig.SeqLengths[0] <= overlapSize*2 { //query in one piece -- fully covered
				seqSet.SetIgnore(contig.Parts[0], true)
			}
			queryStart := contig.Offsets[0]
			queryEnd := queryStart + contig.Lengths[0]
			printLock.Lock()
			for i, part := range contig.Parts[1:] {
				//TODO: what about first seed in this part, find best spot in query? Yes.
				id := i + 1
				rc := "+"
				start := contig.Offsets[id]
				end := start + contig.Lengths[id]
				if contig.ReverseComplement[0] != contig.ReverseComplement[id] {
					rc = "-"
				}
				if contig.SeqLengths[id] <= overlapSize { //this sequence has been fully covered
					seqSet.SetIgnore(part, true)
				}
				ident, _ := contig.Matches[i].GetBasesCovered(k)
				fmt.Printf("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%v\t0\t255\n", seqSet.GetName(contig.Parts[0]), contig.SeqLengths[0], queryStart, queryEnd, rc, seqSet.GetName(part), contig.SeqLengths[id], start, end, ident)
			}
			printLock.Unlock()
		}
	}
	done <- true
}
