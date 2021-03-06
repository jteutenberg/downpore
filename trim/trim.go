package trim

import (
	"fmt"
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/util"
	"log"
	"strings"
	"sync"
)

type Trimmer struct {
	originalFront []sequence.Sequence
	originalBack  []sequence.Sequence

	frontAdapters    []*seeds.SeedSequence
	backAdapters     []*seeds.SeedSequence
	frontAdapterSets []*util.IntSet
	backAdapterSets  []*util.IntSet
	pairsFront       []int
	pairsBack        []int
	index            *seeds.SeedIndex
	k                int
	chunkSize        int

	lock        *sync.Mutex
	frontCounts []int
	backCounts  []int
	noCount     int
	seenCount   int

	midThreshold  int
	extraEdgeTrim int
	extraMidTrim  int
	keepSplits    bool
	tagAdapters   bool
	requirePairs  bool
	verbosity     int
}

type sequenceSplit struct {
	id     int
	aEnd   int
	bStart int
}

//NewTrimmer creates a new Trimmer object that will search for and trim the given adapters.
func NewTrimmer(frontAdapters, backAdapters []sequence.Sequence, k int) *Trimmer {
	var lock sync.Mutex
	t := Trimmer{originalFront: frontAdapters, originalBack: backAdapters, frontAdapters: make([]*seeds.SeedSequence, 0, len(frontAdapters)), backAdapters: make([]*seeds.SeedSequence, 0, len(backAdapters)), frontAdapterSets: make([]*util.IntSet, 0, len(frontAdapters)), backAdapterSets: make([]*util.IntSet, 0, len(backAdapters)), index: nil, k: k, frontCounts: nil, backCounts: nil, noCount: 0, lock: &lock, verbosity: 1}
	(&t).setupIndex()
	t.SetTrimParams(85, 5, 50, 1000, false, true, false)
	return &t
}

func (t *Trimmer) setupIndex() {
	t.frontAdapters = t.frontAdapters[:0]
	t.frontAdapterSets = t.frontAdapterSets[:0]
	t.backAdapters = t.backAdapters[:0]
	t.backAdapterSets = t.backAdapterSets[:0]
	t.index = seeds.NewSeedIndex(uint(t.k))
	//extract all seeds and set up the index
	for _, s := range t.originalFront {
		t.frontAdapters = append(t.frontAdapters, t.index.NewAllSeedSequence(s))
		set := util.NewIntSet()
		kmers := s.ShortKmers(t.k, true)
		t.index.GetSeedsFromKmers(kmers, set)
		t.frontAdapterSets = append(t.frontAdapterSets, set)
	}
	for _, s := range t.originalBack {
		t.backAdapters = append(t.backAdapters, t.index.NewAllSeedSequence(s))
		set := util.NewIntSet()
		kmers := s.ShortKmers(t.k, true)
		t.index.GetSeedsFromKmers(kmers, set)
		t.backAdapterSets = append(t.backAdapterSets, set)
	}
	t.frontCounts = make([]int, len(t.originalFront))
	t.backCounts = make([]int, len(t.originalBack))
	//also prepare the pairs based on adapter names
	pairID := 1
	t.pairsFront = make([]int, len(t.originalFront))
	t.pairsBack = make([]int, len(t.originalBack))
	for i, _ := range t.originalBack {
		t.pairsBack[i] = -1
	}
	for i, a := range t.originalFront {
		name := a.GetName()
		t.pairsFront[i] = -1
		for j, b := range t.originalBack {
			if b.GetName() == name {
				t.pairsFront[i] = pairID
				t.pairsBack[j] = pairID
				pairID++
				break
			}
		}
	}
}

//LoadTrimmer creates a new trimmer loading adapters from fasta/fastq files
func LoadTrimmer(frontAdapters, backAdapters string, k int) *Trimmer {
	fronts := make([]sequence.Sequence, 0, 100)
	backs := make([]sequence.Sequence, 0, 100)
	seqSet := sequence.NewFastaSequenceSet(frontAdapters, 0, 1, false, false)
	seqs := seqSet.GetSequences()
	for seq := range seqs {
		fronts = append(fronts, seq)
	}
	seqSet = sequence.NewFastaSequenceSet(backAdapters, 0, 1, false, false)
	seqs = seqSet.GetSequences()
	for seq := range seqs {
		backs = append(backs, seq)
	}
	return NewTrimmer(fronts, backs, k)
}

//SetVerbosity alters the level of logging output. 0=no output, 1=default, 2=additional output.
func (t *Trimmer) SetVerbosity(level int) {
	t.verbosity = level
}

//SetTrimParams updates the parameters to be used by any following calls to Trim.
func (t *Trimmer) SetTrimParams(midThreshold, extraEdgeTrim, extraMidTrim int, chunkSize int, keepSplits, tagAdapters, requirePairs bool) {
	t.midThreshold = midThreshold
	t.extraEdgeTrim = extraEdgeTrim
	t.extraMidTrim = extraMidTrim
	t.keepSplits = keepSplits
	t.tagAdapters = tagAdapters
	t.requirePairs = requirePairs
	t.chunkSize = chunkSize
}

//Trim searches all sequences for adapters and updates the SequenceSet's offsets accordingly
//Any reads that are split will be disabled in the sequence set and replaced with a new pair (or discarded, depending on current parameters)
func (t *Trimmer) Trim(seqs sequence.SequenceSet, numWorkers int) {
	//read in all the sequences
	ss := seqs.GetSequences()
	//setup workers to build front and back indices
	if t.verbosity > 0 {
		log.Println("Trimming ends and indexing all sequences against", len(t.frontAdapters), "adapters...")
	}
	done := make(chan bool, numWorkers)
	for i := 0; i < numWorkers; i++ {
		go t.trimWorker(seqs, ss, done)
	}
	//wait for completion
	for i := 0; i < numWorkers; i++ {
		<-done
	}
	//add some sequences to the index
	ss = seqs.GetSequences()
	longestAdapter := 100 //bases in longest adapter, with a bit of padding
	edgeSize := 150       //bases to search for early and late adapters
	minSeeds := 4         //minimum number of seeds required to make an splitting match with an adapter
	totalCount := 0
	totalBases := 0

	splits := make([]*sequenceSplit, seqs.Size()+1)
	ids := make([]int, 0, 100)
	var maxID int
	for seq := range ss {
		//put the remaining centre of the sequence into the index for later querying
		totalBases += seq.Len()-edgeSize*2
		for i := edgeSize; i < seq.Len()-edgeSize-longestAdapter; i += t.chunkSize - longestAdapter { //100 is minimum non-edge sequence size, and max expected adapter size
			if i > seq.Len()-(t.chunkSize*3)/2-edgeSize {
				//add the entire remainder
				seedSeq := t.index.NewSeedSequence(seq.SubSequence(i, seq.Len()-edgeSize))
				totalCount += seedSeq.GetNumSeeds()
				t.index.AddSequence(seedSeq)
				break
			} else {
				//just a chunk
				endPoint := i + t.chunkSize
				if endPoint >= seq.Len()-edgeSize {
					endPoint = seq.Len() - edgeSize
				}
				seedSeq := t.index.NewSeedSequence(seq.SubSequence(i, endPoint))
				totalCount += seedSeq.GetNumSeeds()
				if seedSeq.GetNumSeeds() >= minSeeds {
					t.index.AddSequence(seedSeq)
				}
			}
		}
		//refresh the index every few hundred thousand sequences? Or just count #seeds
		if totalCount > 300000000 { //300M
			t.index.IndexSequences(numWorkers)
			//now query for internal matches
			if t.verbosity > 0 {
				log.Println("Searching", int(totalBases/1000000),"MB of sequences for splitting based on", len(t.frontAdapters), "adapters")
			}

			for i, ad := range t.frontAdapters {
				go t.findSplit(ad, t.frontAdapterSets[i], splits, &ids, &maxID, seqs, done)
			}
			for i := 0; i < len(t.frontAdapters); i++ {
				<-done
			}
			totalCount = 0
			totalBases = 0
			//clear the sequences
			t.setupIndex()
		}
	}
	//and do the final splits too
	if totalCount > 0 {
		t.index.IndexSequences(numWorkers)
		if t.verbosity > 0 {
			log.Println("Searching", int(totalBases/1000000),"MB of sequences for splitting based on", len(t.frontAdapters), "adapters")
		}
		for i, ad := range t.frontAdapters {
			go t.findSplit(ad, t.frontAdapterSets[i], splits, &ids, &maxID, seqs, done)
		}
		for i := 0; i < len(t.frontAdapters); i++ {
			<-done
		}
	}
	if t.verbosity > 0 {
		log.Println(len(ids), "sequences require splitting")
	}
	//get the full (now edge-trimmed) sequences that need splitting
	ss = seqs.GetSequencesByID(ids)
	splitSeqs := make([]sequence.Sequence, maxID+1, maxID+1)
	for s := range ss {
		splitSeqs[s.GetID()] = s
	}
	for _, id := range ids {
		split := splits[id]
		if split == nil {
			continue
		}
		seq := splitSeqs[id]
		if t.keepSplits {
			report := fmt.Sprint("Splitting read ", split.id, " into")
			//test for invalid splits. This occurs when a second internal adapter reduces one side to less than minimum sequence length
			if split.aEnd > edgeSize {
				seqs.AddSequence(seq.SubSequence(0, split.aEnd), seqs.GetName(split.id)+"_(left)")
				report += fmt.Sprint(": 0 - ",split.aEnd," and ")
			} else {
				report += " ignored short left hand side and "
			}
			if seq.Len()-split.bStart > edgeSize {
				seqs.AddSequence(seq.SubSequence(split.bStart, seq.Len()), seqs.GetName(split.id)+"_(right)")
				report += fmt.Sprint(split.bStart," - ",seq.Len())
			} else {
				report += " ignored short right hand side"
			}
			if t.verbosity > 1 {
				log.Println(report)
				if split.aEnd >= 0 && split.bStart < seq.Len() && split.bStart-split.aEnd-t.extraMidTrim*2 <= longestAdapter {
					log.Println(seq.SubSequence(split.aEnd+t.extraMidTrim, split.bStart-t.extraMidTrim).String())
				}
			}
		}
		seqs.SetIgnore(split.id, true)
	}
}

//PrintStats outputs to stderr all adapters present (reduced after DetermineAdapters is called) and the percentage of reads containing them (modified after a call to Trim)
func (t *Trimmer) PrintStats(seqs sequence.SequenceSet) {
	for i, count := range t.frontCounts {
		log.Println("Front adapter:", t.originalFront[i].GetName(), "\t", (count*100)/t.seenCount, "%")
	}
	for i, count := range t.backCounts {
		log.Println("Back adapter:", t.originalBack[i].GetName(), "\t", (count*100)/t.seenCount, "%")
	}
	log.Println((t.noCount*100)/t.seenCount, "% with no adapters found.")
}

//DetermineAdapters runs through the first numReads reads, discarding adapters from the trimmer's set that
//do not obtain any good matches
func (t *Trimmer) DetermineAdapters(seqs sequence.SequenceSet, numReads int, threshold int, numWorkers int) {
	ss := seqs.GetNSequencesFrom(0, numReads)
	frontEnabled := make([]bool, len(t.frontAdapters), len(t.frontAdapters))
	backEnabled := make([]bool, len(t.backAdapters), len(t.backAdapters))
	//start workers
	done := make(chan bool, numWorkers)
	for i := 0; i < numWorkers; i++ {
		go t.checkAdapterWorker(seqs, frontEnabled, backEnabled, threshold, ss, done)
	}
	//wait for completion
	for i := 0; i < numWorkers; i++ {
		<-done
	}
	count := 0
	for _, en := range frontEnabled {
		if en {
			count++
		}
	}
	if t.verbosity > 0 {
		log.Println(count, "/", len(frontEnabled), "front adapters identified with high identity matches.")
	}
	for i := len(frontEnabled) - 1; i >= 0; i-- {
		if frontEnabled[i] {
			if t.verbosity > 0 {
				log.Println(" -", t.originalFront[i].GetName())
			}
		} else {
			t.originalFront[i] = t.originalFront[len(t.originalFront)-1]
			t.originalFront = t.originalFront[:len(t.originalFront)-1]
		}
	}
	count = 0
	for _, en := range backEnabled {
		if en {
			count++
		}
	}
	if t.verbosity > 0 {
		log.Println(count, "/", len(backEnabled), "back adapters identified with high identity matches.")
	}
	for i := len(backEnabled) - 1; i >= 0; i-- {
		if backEnabled[i] {
			if t.verbosity > 0 {
				log.Println(" -", t.originalBack[i].GetName())
			}
		} else {
			t.originalBack[i] = t.originalBack[len(t.originalBack)-1]
			t.originalBack = t.originalBack[:len(t.originalBack)-1]
		}
	}
	t.setupIndex()
}

func (t *Trimmer) isNewFullMatch(kmerSet *util.IntSet, seq sequence.Sequence, threshold int, adapters []*seeds.SeedSequence, adapterSets []*util.IntSet, enabled []bool) {
	var seedSeq *seeds.SeedSequence
	for i, adapter := range adapterSets {
		if enabled[i] {
			continue //we already know this is a good adapter
		}
		hits := kmerSet.CountIntersection(adapter)
		minHits := adapter.Size() / 2 //50% shared kmers minimum (test actual identity later)
		if hits >= minHits {
			//test their ordering
			if seedSeq == nil {
				seedSeq = t.index.NewSeedSequence(seq)
			}
			ms := seedSeq.Match(adapters[i], adapter, kmerSet, int(minHits-1), t.k) //require the same high number of in-order matches
			if ms != nil {
				for _, m := range ms {
					if len(m.MatchA) >= int(minHits) {
						identity, _ := m.GetBasesCovered(t.k)
						if (identity*100)/adapters[i].Len() >= threshold {
							enabled[i] = true
						}
					}
				}
			}
		}
	}
}

func (t *Trimmer) findMatches(kmerSet *util.IntSet, seq sequence.Sequence, adapters []*seeds.SeedSequence, adapterSets []*util.IntSet, counts []int) (int, int, bool, int) {
	var seedSeq *seeds.SeedSequence
	earliest := seq.Len()
	latest := 0
	found := false
	var bestMatch int
	var bestIdent int
	var barcoded bool
	var ambiguous bool
	for i, adapter := range adapterSets {
		hits := kmerSet.CountIntersection(adapter)
		fraction := (hits * 10) / adapter.Size()
		if fraction >= 2 || hits >= 3 { //enough matching to go to next test
			if seedSeq == nil {
				seedSeq = t.index.NewSeedSequence(seq)
			}
			//have enough k-mers, check their ordering
			ms := seedSeq.Match(adapters[i], adapter, kmerSet, 3, t.k)
			if ms != nil {
				for _, m := range ms {
					if len(m.MatchA) >= 3 {
						identity, _ := m.GetBasesCovered(t.k)
						identity = (identity * 100) / adapters[i].Len()
						isBarcode := strings.HasPrefix(adapters[i].GetName(), "Barcode")
						if !barcoded && isBarcode {
							barcoded = true
							bestIdent = identity
							bestMatch = i
						} else if barcoded {
							if isBarcode {
								//test for ambiguity
								delta := identity - bestIdent
								ambiguous = delta < 5 && delta > -5
								if identity > bestIdent {
									bestIdent = identity
									bestMatch = i
								}
							}
						} else if identity > bestIdent {
							bestIdent = identity
							bestMatch = i
						}

						start := seedSeq.GetSeedOffset(m.MatchB[0], t.k) + adapters[i].GetSeedOffset(m.MatchA[0], t.k)
						end := seedSeq.GetSeedOffset(m.MatchB[len(m.MatchB)-1], t.k) + adapters[i].GetSeedOffsetFromEnd(m.MatchA[len(m.MatchA)-1], t.k)
						if start < earliest {
							if start < 0 {
								start = 0
							}
							earliest = start
						}
						if end > latest {
							if end > seq.Len() {
								end = seq.Len()
							}
							latest = end
						}
						found = true
						t.lock.Lock()
						if i > len(counts) {
							log.Println("Over count! ", i, " adapter out of ", len(counts))
						}
						counts[i]++
						t.lock.Unlock()
					}
				}
			}
		}
	}
	if ambiguous {
		//trim, but pretend we didn't see an adapter
		return earliest, latest, false, 0
	}
	return earliest, latest, found, bestMatch
}

func (t *Trimmer) checkAdapterWorker(set sequence.SequenceSet, frontEnabled, backEnabled []bool, threshold int, seqs <-chan sequence.Sequence, done chan<- bool) {
	kmerSet := util.NewIntSet()
	edgeSize := 150 //bases to search for early and late adapters
	for seq := range seqs {
		if seq.Len() < edgeSize+50 {
			continue
		}
		frontSeq := seq.SubSequence(0, edgeSize)
		backSeq := seq.SubSequence(seq.Len()-edgeSize, seq.Len())
		kmers := frontSeq.ShortKmers(t.k, true)
		kmerSet.Clear()
		t.index.GetSeedsFromKmers(kmers, kmerSet)
		t.isNewFullMatch(kmerSet, frontSeq, threshold, t.frontAdapters, t.frontAdapterSets, frontEnabled)
		kmerSet.Clear()
		kmers = backSeq.ShortKmers(t.k, true)
		t.index.GetSeedsFromKmers(kmers, kmerSet)
		t.isNewFullMatch(kmerSet, backSeq, threshold, t.backAdapters, t.backAdapterSets, backEnabled)
	}
	done <- true
}

func (t *Trimmer) trimWorker(set sequence.SequenceSet, seqs <-chan sequence.Sequence, done chan<- bool) {
	kmerSet := util.NewIntSet()
	edgeSize := 150 //bases to search for early and late adapters. TODO: pass in as an argument, ensure consistency with the internal splitting
	for seq := range seqs {
		if seq.Len() < edgeSize+50 {
			continue
		}
		frontSeq := seq.SubSequence(0, edgeSize)
		backSeq := seq.SubSequence(seq.Len()-edgeSize, seq.Len())
		kmers := frontSeq.ShortKmers(t.k, true) //more sensitive
		//manually check first 150 vs front adapters and last 150 vs back adapters
		kmerSet.Clear()
		t.index.GetSeedsFromKmers(kmers, kmerSet)
		_, start, foundStart, matchIndex := t.findMatches(kmerSet, frontSeq, t.frontAdapters, t.frontAdapterSets, t.frontCounts)
		kmerSet.Clear()
		kmers = backSeq.ShortKmers(t.k, true) //more sensitive
		t.index.GetSeedsFromKmers(kmers, kmerSet)
		end, _, foundEnd, backMatchIndex := t.findMatches(kmerSet, backSeq, t.backAdapters, t.backAdapterSets, t.backCounts)

		//if pairs are required, test whether these exist
		if t.requirePairs {
			f := -1
			b := -1
			if foundStart {
				f = t.pairsFront[matchIndex]
			}
			if foundEnd {
				b = t.pairsBack[backMatchIndex]
			}
			if f != b { //one of them wants to pair, but the other doesn't match
				//disable the adapter match, but still trim. Same as an ambiguous adapter match.
				foundStart = false
				foundEnd = false
			}
		}

		t.lock.Lock()
		t.seenCount++
		if !foundStart {
			t.noCount++
		}
		t.lock.Unlock()
		start += t.extraEdgeTrim
		end = edgeSize - end + t.extraEdgeTrim //convert to trim amount
		if start+end+10 >= seq.Len() {
			set.SetIgnore(seq.GetID(), true)
		} else {
			if foundStart {
				set.SetFrontTrim(seq.GetID(), start)
				if t.tagAdapters {
					set.SetName(seq.GetID(), t.frontAdapters[matchIndex].GetName()+"_"+set.GetName(seq.GetID()))
				}
			} else if end > start && start > 0 {
				//trim off ambiguous adapters too
				set.SetFrontTrim(seq.GetID(), start)
			}
			if foundEnd || (end > start && end < seq.Len()) {
				set.SetBackTrim(seq.GetID(), end)
			}
		}
	}
	done <- true
}

func (t *Trimmer) findSplit(ad *seeds.SeedSequence, adSet *util.IntSet, splits []*sequenceSplit, ids *[]int, maxID *int, seqs sequence.SequenceSet, done chan<- bool) {
	//edgeSize := 150
	minSeqLength := 500 //splits need to be longer than this

	minMatch := ad.GetNumSeeds() / 5
	ms := t.index.Matches(ad, 0.2)
	for _, index := range ms {
		target := t.index.GetSeedSequence(uint(index))
		targetSet := t.index.GetSeedSet(uint(index))
		matches := target.Match(ad, adSet, targetSet, minMatch, t.k)
		if matches != nil {
			for _, match := range matches {
				identity, _ := match.GetBasesCovered(t.k)
				if (identity*100)/ad.Len() < t.midThreshold {
					continue
				}
				t.lock.Lock()

				id := target.GetID()
				if id < 0 || id >= len(splits) {
					log.Println("Warning: unexpected sequence for splitting, id: ", id, "/", len(splits))
					continue
				}
				frontTrim := seqs.GetFrontTrim(id)
				backTrim := seqs.GetBackTrim(id)

				start := target.GetOffset() + target.GetSeedOffset(match.MatchB[0], t.k) - ad.GetSeedOffset(match.MatchA[0], t.k)
				seqLen := target.GetOffset() + target.Len() + target.GetInset() - backTrim
				if start < minSeqLength+frontTrim { //just crop the front off
					newTrim := start+ad.Len()+t.extraMidTrim
					if newTrim+minSeqLength < seqLen {
						if newTrim > frontTrim {
							seqs.SetFrontTrim(target.GetID(), newTrim)
							if splits[id] != nil { //update the existing split
								splits[id].aEnd -= (newTrim-frontTrim)
								splits[id].bStart -= (newTrim-frontTrim)
							}
						}
						if t.tagAdapters {
							seqs.SetName(target.GetID(), ad.GetName()+"_"+seqs.GetName(target.GetID()))
						}
					} else {
						splits[id] = nil //in case of existing split that is no longer valid
						seqs.SetIgnore(target.GetID(), true)
					}
				} else if start+minSeqLength+ad.Len() > seqLen { //crop off the tail
					newTrim := seqLen-start+t.extraMidTrim
					if newTrim > backTrim {
						seqs.SetBackTrim(target.GetID(), newTrim)
					}
				} else {
					//prepare the split. Compensate for the trim that will be applied.
					//this is because we are working with the seed sequences that are still in a pre-trim state, but will apply the split to the trimmed sequence
					if splits[id] != nil {
						if splits[id].aEnd > start-t.extraMidTrim-frontTrim {
							splits[id].aEnd = start - t.extraMidTrim - frontTrim
						}
						if splits[id].bStart < start+ad.Len()+t.extraMidTrim-frontTrim {
							splits[id].bStart = start + ad.Len() + t.extraMidTrim - frontTrim
						}
					} else {
						if t.verbosity > 2 {
							log.Println(match.LongString(t.index))
						}
						splits[id] = &sequenceSplit{id: id, aEnd: start - t.extraMidTrim - frontTrim, bStart: start + ad.Len() + t.extraMidTrim - frontTrim}
						*ids = append(*ids, id)
						if id > *maxID {
							*maxID = id
						}
					}
				}
				t.lock.Unlock()
			}
		}
	}
	done <- true
}
