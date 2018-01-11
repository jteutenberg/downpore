package trim

import (
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
	index            *seeds.SeedIndex
	k                int

	lock        *sync.Mutex
	frontCounts []int
	backCounts  []int
	noCount     int

	endThreshold  int
	midThreshold  int
	absThreshold  int
	extraEdgeTrim int
	extraMidTrim  int
	keepSplits    bool
	tagAdapters   bool
	verbosity     int
}

type sequenceSplit struct {
	id     int
	aEnd   int
	bStart int
}

func NewTrimmer(frontAdapters, backAdapters []sequence.Sequence, k int) *Trimmer {
	var lock sync.Mutex
	t := Trimmer{originalFront: frontAdapters, originalBack: backAdapters, frontAdapters: make([]*seeds.SeedSequence, 0, len(frontAdapters)), backAdapters: make([]*seeds.SeedSequence, 0, len(backAdapters)), frontAdapterSets: make([]*util.IntSet, 0, len(frontAdapters)), backAdapterSets: make([]*util.IntSet, 0, len(backAdapters)), index: nil, k: k, frontCounts: nil, backCounts: nil, noCount: 0, lock: &lock, verbosity: 1}
	(&t).setupIndex()
	t.SetTrimParams(75, 75, 25, 5, 50, false, true)
	return &t
}

func (t *Trimmer) setupIndex() {
	t.index = seeds.NewSeedIndex(uint(t.k))
	//extract all seeds and set up the index
	for _, s := range t.originalFront {
		t.frontAdapters = append(t.frontAdapters, t.index.NewAllSeedSequence(s))
		set := util.NewIntSet()
		t.frontAdapterSets = append(t.frontAdapterSets, set)
		kmers := s.ShortKmers(t.k)
		for _, kmer := range kmers {
			set.Add(uint(kmer))
		}

	}
	for _, s := range t.originalBack {
		t.backAdapters = append(t.backAdapters, t.index.NewAllSeedSequence(s))
		set := util.NewIntSet()
		t.backAdapterSets = append(t.backAdapterSets, set)
		kmers := s.ShortKmers(t.k)
		for _, kmer := range kmers {
			set.Add(uint(kmer))
		}
	}
	t.frontCounts = make([]int, len(t.originalFront), len(t.originalFront))
	t.backCounts = make([]int, len(t.originalBack), len(t.originalBack))
}

//LoadTrimmer creates a new trimmer loading adapters from fasta/fastq files
func LoadTrimmer(frontAdapters, backAdapters string, k int) *Trimmer {
	fronts := make([]sequence.Sequence, 0, 100)
	backs := make([]sequence.Sequence, 0, 100)
	seqSet := sequence.NewFastaSequenceSet(frontAdapters, 0)
	seqs := seqSet.GetSequences()
	for seq := range seqs {
		fronts = append(fronts, seq)
	}
	seqSet = sequence.NewFastaSequenceSet(backAdapters, 0)
	seqs = seqSet.GetSequences()
	for seq := range seqs {
		backs = append(backs, seq)
	}
	return NewTrimmer(fronts, backs, k)
}

func (t *Trimmer) SetVerbosity(level int) {
	t.verbosity = level
}

func (t *Trimmer) SetTrimParams(endThreshold, midThreshold, absThreshold, extraEdgeTrim, extraMidTrim int, keepSplits bool, tagAdapters bool) {
	t.endThreshold = endThreshold
	t.midThreshold = midThreshold
	t.absThreshold = absThreshold
	t.extraEdgeTrim = extraEdgeTrim
	t.extraMidTrim = extraMidTrim
	t.keepSplits = keepSplits
	t.tagAdapters = tagAdapters
}

//Trim searches all sequences for adapters and updates the SequenceSet's offsets accordingly
//Any reads that are split will be disabled in the sequence set and replaced with a new pair
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

	edgeSize := 150
	minSeqLength := 1000 //splits need to be longer than this

	//now query for internal matches
	if t.verbosity > 0 {
		log.Println("Searching", t.index.GetNumSequences(), "sequences for splitting based on", len(t.frontAdapters), "adapters")
	}
	splits := make([]*sequenceSplit, t.index.GetNumSequences(), t.index.GetNumSequences())
	ids := make([]int, 0, 100)
	var maxID int
	for i, ad := range t.frontAdapters {
		minMatch := ad.GetNumSeeds() / 3
		ms := t.index.Matches(ad, 0.5)
		//fire off workers to handle the matches
		count := 0
		for _, index := range ms {
			target := t.index.GetSeedSequence(uint(index))
			match := target.Match(ad, t.frontAdapterSets[i], minMatch, t.k)
			if match != nil && len(match.MatchA) >= minMatch {
				identity, _ := match.GetBasesCovered(t.k)
				if identity < t.absThreshold && (identity*100)/ad.Len() < t.midThreshold {
					continue
				}
				start := target.GetOffset() + target.GetSeedOffset(match.MatchB[0], t.k) - ad.GetSeedOffset(match.MatchA[0], t.k)
				seqLen := target.GetOffset() + target.Len() + target.GetInset()
				if start < minSeqLength { //just crop the front off
					if start+ad.Len()+edgeSize < seqLen {
						seqs.SetFrontTrim(target.GetID(), start+ad.Len()+t.extraMidTrim)
						if t.tagAdapters {
							seqs.SetName(target.GetID(), ad.GetName()+" "+seqs.GetName(target.GetID()))
						}
					} else {
						seqs.SetIgnore(target.GetID(), true)
					}
				} else if start+minSeqLength+ad.Len() > seqLen { //crop off the tail
					seqs.SetBackTrim(target.GetID(), seqLen-start+t.extraMidTrim)
				} else {
					//prepare the split. Compensate for the trim that will be applied.
					id := target.GetID()
					futureTrim := seqs.GetFrontTrim(id)
					if splits[id] != nil {
						if splits[id].aEnd > start-t.extraMidTrim-futureTrim {
							splits[id].aEnd = start - t.extraMidTrim - futureTrim
						}
						if splits[id].bStart < start+ad.Len()+t.extraMidTrim-futureTrim {
							splits[id].bStart = start + ad.Len() + t.extraMidTrim - futureTrim
						}
					} else {
						splits[id] = &sequenceSplit{id: id, aEnd: start - t.extraMidTrim - futureTrim, bStart: start + ad.Len() + t.extraMidTrim - futureTrim}
						ids = append(ids, id)
						if id > maxID {
							maxID = id
						}
					}
				}
				count++
			}
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
		seq := splitSeqs[id]
		//Note: if it is 3+ combined sequences it only keeps the first and last
		if t.keepSplits {
			if t.verbosity > 1 {
				log.Println("Splitting read ", split.id, " into: 0 -", split.aEnd, "and", split.bStart, "-", seq.Len())
				log.Println(seq.SubSequence(split.aEnd+edgeSize, split.bStart-edgeSize).String())
			}
			seqs.AddSequence(seq.SubSequence(0, split.aEnd), seqs.GetName(split.id)+" (left)")
			seqs.AddSequence(seq.SubSequence(split.bStart, seq.Len()), seqs.GetName(split.id)+" (right)")
		}
		seqs.SetIgnore(split.id, true)
	}
}

func (t *Trimmer) PrintStats(seqs sequence.SequenceSet) {
	for i, count := range t.frontCounts {
		log.Println("Front adapter:", t.originalFront[i].GetName(), "\t", (count*100)/seqs.Size(), "%")
	}
	for i, count := range t.backCounts {
		log.Println("Back adapter:", t.originalBack[i].GetName(), "\t", (count*100)/seqs.Size(), "%")
	}
	log.Println((t.noCount*100)/seqs.Size(), "% with no adapters found.")
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
	t.frontAdapters = t.frontAdapters[:0]
	t.frontAdapterSets = t.frontAdapterSets[:0]
	t.backAdapters = t.backAdapters[:0]
	t.backAdapterSets = t.backAdapterSets[:0]
	t.setupIndex()
}

func (t *Trimmer) isNewFullMatch(kmerSet *util.IntSet, seq sequence.Sequence, threshold int, adapters []*seeds.SeedSequence, adapterSets []*util.IntSet, enabled []bool) {
	var seedSeq *seeds.SeedSequence
	for i, adapter := range adapterSets {
		if enabled[i] {
			continue //we already know this is a good adapter
		}
		hits := kmerSet.CountIntersection(adapter)
		minHits := (adapter.Size() * 7) / 10 //70% shared kmers minimum
		if hits >= minHits {
			//test their ordering
			if seedSeq == nil {
				seedSeq := t.index.NewSeedSequence(seq, 0, nil)
				m := seedSeq.Match(adapters[i], adapter, int(minHits-1), t.k) //require the same high number of in-order matches
				if m != nil && len(m.MatchA) >= int(minHits) {
					identity, _ := m.GetBasesCovered(t.k)
					if (identity*100)/adapters[i].Len() >= threshold {
						enabled[i] = true
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
	var bestHits int
	var barcoded bool
	for i, adapter := range adapterSets {
		hits := kmerSet.CountIntersection(adapter)
		fraction := (hits * 10) / adapter.Size()
		if fraction >= 3 || hits > 4 { //enough matching to go to next test
			if seedSeq == nil {
				seedSeq = t.index.NewSeedSequence(seq, 0, nil)
			}
			//have enough k-mers, check their ordering
			m := seedSeq.Match(adapters[i], adapter, 4, t.k)
			if m != nil && len(m.MatchA) >= 4 {
				if !barcoded {
					name := adapters[i].GetName()
					if strings.HasPrefix(name, "Barcode") {
						barcoded = true
						bestMatch = i
					} else if len(m.MatchA) > bestHits {
						bestHits = len(m.MatchA)
						bestMatch = i
					}
				}
				//loose criteria is fine. We want high recall, precision isn't so important.
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
				counts[i]++
				t.lock.Unlock()
			}
		}
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
		kmers := seq.ShortKmers(t.k)
		kmerSet.Clear()
		for _, k := range kmers[:edgeSize] {
			kmerSet.Add(uint(k))
		}
		t.isNewFullMatch(kmerSet, seq.SubSequence(0, edgeSize), threshold, t.frontAdapters, t.frontAdapterSets, frontEnabled)
		kmerSet.Clear()
		for _, k := range kmers[len(kmers)-edgeSize:] {
			kmerSet.Add(uint(k))
		}
		t.isNewFullMatch(kmerSet, seq.SubSequence(seq.Len()-edgeSize, seq.Len()), threshold, t.backAdapters, t.backAdapterSets, backEnabled)
	}
	done <- true
}

func (t *Trimmer) trimWorker(set sequence.SequenceSet, seqs <-chan sequence.Sequence, done chan<- bool) {
	kmerSet := util.NewIntSet()
	edgeSize := 150       //bases to search for early and late adapters
	longestAdapter := 100 //bases in longest adapter, with a bit of padding
	minSeeds := 4         //minimum number of seeds required to make an splitting match with an adapter
	for seq := range seqs {
		if seq.Len() < edgeSize+50 {
			continue
		}
		kmers := seq.ShortKmers(t.k) //more sensitive
		//manually check first 150 vs front adapters and last 150 vs back adapters
		kmerSet.Clear()
		for _, k := range kmers[:edgeSize] {
			kmerSet.Add(uint(k))
		}
		_, start, foundStart, matchIndex := t.findMatches(kmerSet, seq.SubSequence(0, edgeSize), t.frontAdapters, t.frontAdapterSets, t.frontCounts)
		kmerSet.Clear()
		for _, k := range kmers[len(kmers)-edgeSize:] {
			kmerSet.Add(uint(k))
		}
		end, _, foundEnd, _ := t.findMatches(kmerSet, seq.SubSequence(seq.Len()-edgeSize, seq.Len()), t.backAdapters, t.backAdapterSets, t.backCounts)
		if !foundStart {
			t.lock.Lock()
			t.noCount++
			t.lock.Unlock()
		}
		start += t.extraEdgeTrim
		end = edgeSize - end + t.extraEdgeTrim //convert to trim amount
		if start+end+10 >= seq.Len() {
			set.SetIgnore(seq.GetID(), true)
		} else {
			if foundStart {
				set.SetFrontTrim(seq.GetID(), start)
				if t.tagAdapters {
					set.SetName(seq.GetID(), t.frontAdapters[matchIndex].GetName()+" "+set.GetName(seq.GetID()))
				}
			}
			if foundEnd {
				set.SetBackTrim(seq.GetID(), end)
			}
			//put the remaining centre of the sequence into the index for later querying
			chunkSize := 1000
			for i := edgeSize; i < seq.Len()-edgeSize-longestAdapter; i += chunkSize - longestAdapter { //100 is minimum non-edge sequence size, and max expected adapter size
				if i > seq.Len()-(chunkSize*3)/2-edgeSize {
					//add the entire remainder
					seedSeq := t.index.NewSeedSequence(seq.SubSequence(i, seq.Len()-edgeSize), 0, nil)
					t.index.AddSequence(seedSeq)
					break
				} else {
					//just a chunk
					endPoint := i + chunkSize
					if endPoint >= seq.Len()-edgeSize {
						endPoint = seq.Len() - edgeSize
					}
					seedSeq := t.index.NewSeedSequence(seq.SubSequence(i, endPoint), 0, nil)
					if seedSeq.GetNumSeeds() >= minSeeds {
						t.index.AddSequence(seedSeq)
					}
				}
			}
		}
	}
	done <- true
}
