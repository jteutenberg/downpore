package commands

import (
	"bufio"
	"fmt"
	"github.com/jteutenberg/downpore/sequence"
	"github.com/jteutenberg/downpore/seeds"
	"github.com/jteutenberg/downpore/util/formats"
	"github.com/jteutenberg/downpore/util/sequtil"
	"math"
	"os"
	"runtime"
	"runtime/debug"
	"sort"
	"strings"
	"sync"
)

type kmersCommand struct {
	args  map[string]string
	alias map[string]string
	desc  map[string]string

	lock *sync.Mutex
	goodQuality []int32
	badQuality []int32
	kmerGoodCounts []int32 //number of correct instances of this kmer
	kmerBadCounts []int32
	kmerGoodQuality []int32 //quality over all correct instances
	kmerBadQuality []int32
}

func NewKmersCommand() Command {
	args, alias, desc := MakeArgs(
		[]string{"input", "alignment","reference", "training","training_alignment","training_ref","k","map_size","num_workers"},
		[]string{"", "", "","","","","10","100","4"},
		[]string{"Reads input file", "SAM input file", "Reference fasta file","Training input file","SAM training file","Training reference fasta file","K-mer size","Dimensions for heatmaps","Number of worker threads to use"})
	var lock sync.Mutex
	return &kmersCommand{args: args, alias: alias, desc: desc, lock: &lock}
}
func (com *kmersCommand) GetName() string {
	return "kmers"
}
func (com *kmersCommand) GetArgs() (map[string]string, map[string]string, map[string]string) {
	return com.args, com.alias, com.desc
}

//reurns the new indices, shrunk to remove invalid data
func (com *kmersCommand) getCounts(inputFile, alignmentsFile string, reference string, k, numWorkers int, gCount, bCount, gQuality, bQuality []int32, indices []int) ([]float64,[]int) {
	seqSet := sequence.NewFastaSequenceSet(inputFile, 0, numWorkers, true, false)
	ids := make(map[string]int) //map names to sequence ids
	allSeqs := seqSet.GetSequences()
	nextID := 0
	for s := range allSeqs {
		ids[s.GetName()] = nextID
		nextID++
	}
	alignments := formats.LoadSAM(alignmentsFile)
	prevSeq := ""
	alignIn := make(chan *formats.SAMAlignment, numWorkers*2)
	done := make(chan bool, numWorkers)
	for i := 0; i < numWorkers; i++ {
		go com.alignmentWorker(alignIn, seqSet, ids, reference, k, com.kmerGoodCounts, com.kmerBadCounts, com.kmerGoodQuality, com.kmerBadQuality, done)
	}
	for a := range alignments {
		if a.NameA == prevSeq {
			continue
		}
		prevSeq = a.NameA
		alignIn <- a
	}
	close(alignIn)
	for i := 0; i < numWorkers; i++ {
		<-done
	}

	//convert to ranks, removing irrelevant bits
	values := make([]float64, len(gCount))
	for i, index := range indices {
		if index == math.MaxInt64 {
			continue
		}
		gc := gCount[index]
		bc := bCount[index]
		if gc+bc > 2 { //must occur in the data at least 3 times
			values[i] = float64(gc)/float64(gc+bc)
		} else {
			//remove the index
			indices[i] = math.MaxInt64
		}

	}
	//values and indices are both full length, with indices containing to-remove values
	values, indices = rankify(values,indices)
	return values, indices
}

type datum struct {
	goodCount int16
	badCount int16
	allQ int32
}

//returns: accuracy, quality, rcRatio, indices
func getLongCounts(inputFile, alignmentsFile string, reference string, k int, whitelist map[int]*datum) map[int]*datum {
	seqSet := sequence.NewFastaSequenceSet(inputFile, 0, 4, false, false)
	if whitelist == nil {
		whitelist := make(map[int]*datum,10000000)
		kmerCounts := sequtil.LongKmerOccurrences(seqSet.GetSequences(), k)
		for i, count := range kmerCounts {
			if count > 2 {
				whitelist[i] = nil
			}
		}
		kmerCounts = nil
		runtime.GC()
		debug.FreeOSMemory()
		fmt.Println("Got whitelist from sequences.")
	}

	ids := make(map[string]int) //map names to sequence ids
	allSeqs := seqSet.GetSequences()
	nextID := 0
	for s := range allSeqs {
		ids[s.GetName()] = nextID
		nextID++
	}
	alignments := formats.LoadSAM(alignmentsFile)
	data := make(map[int]*datum, 10000000)
	prevSeq := ""
	for a := range alignments {
		if a.NameA == prevSeq {
			continue
		}
		prevSeq = a.NameA
		seq := <-seqSet.GetNSequencesFrom(ids[a.NameA],1)
		original := seq.String()
		if a.ReverseComplement {
			seq = seq.ReverseComplement()
		}
		s := seq.String()
		q := seq.Quality()

		matches := a.AlignmentString.GetKmerMatches(k)
		prevSPos := 0
		for seqIndex, ok := <-matches; ok; seqIndex, ok = <-matches {
			refIndex, _ := <-matches
			refIndex += a.StartB

			if prevSPos == 0 { //ignore the start and end of the cigar
				prevSPos = seqIndex
			}
			//The match
			sKmer := sequence.KmerValue(original[len(original)-k-seqIndex:len(original)-seqIndex]) //same position in the original
			skip := false
			if whitelist != nil {
				if _, exists := whitelist[sKmer]; !exists {
					skip = true
				}
			}
			var d *datum
			var ok bool
			if !skip {
				if d,ok = data[sKmer]; !ok {
					d = &datum{} //all zero
					data[sKmer] = d
				}
				if reference[refIndex:refIndex+k] == s[seqIndex:seqIndex+k] {
					d.goodCount++
				} else {
					d.badCount++
				}
				if q != nil {
					d.allQ += int32(q[seqIndex+k/2])
				}
			}
			//non-matches up to here
			for prevSPos < seqIndex {
				sKmer = sequence.KmerValue(original[len(s)-k-prevSPos:len(s)-prevSPos])
				skip = false
				if whitelist != nil {
					if _, exists := whitelist[sKmer]; !exists {
						skip = true
					}
				}
				if !skip {
					if d,ok = data[sKmer]; !ok {
						d = &datum{} //all zero
						data[sKmer] = d
					}
					d.badCount++
					if q != nil {
						d.allQ += int32(q[prevSPos+k/2])
					}
				}
				prevSPos++
			}
			prevSPos = seqIndex+1
		}
	}
	fmt.Println("Total k-mers found: ",len(data))
	//clear out some irrelevant data points
	for kmer, d := range data {
		sum := d.goodCount + d.badCount
		if sum <= 2 {
			delete(data,kmer)
		}
	}
	fmt.Println("After removing low frequency: ",len(data))

	return data
}

func getLongCorrelations(data,trainingData map[int]*datum, alignmentsFile string, k int) {

	accuracies := make([]float64, len(data))
	qualities := make([]float64, len(data))
	rcRatios := make([]float64, len(data))
	lex := make([]float64, len(data))
	trained := make([]float64, len(data))
	indices := make([]int, len(data))

	i := 0
	for kmer, d := range data {
		sum := d.goodCount + d.badCount
		if sum > 2 {
			indices[i] = i
			lex[i] = float64(kmer)
			accuracies[i] = float64(d.goodCount)/float64(sum)
			qualities[i] = float64(d.allQ)/float64(sum)
			rc := int(seeds.ReverseComplement(uint(kmer),uint(k)))
			if rcDatum, ok := data[rc]; ok {
				rcRatio := float64(sum) / float64(sum + rcDatum.goodCount + rcDatum.badCount)
				rcRatio= 0.5 - rcRatio
				if rcRatio < 0.0 {
					rcRatio = -rcRatio
				}
				rcRatios[i] = 0.5 - rcRatio
			} else {
				rcRatios[i] = 0.0
			}
			if td, ok := trainingData[kmer]; ok {
				sum := td.goodCount+td.badCount
				if sum > 2 {
					trained[i] = float64(td.goodCount)/float64(sum)
				}
			}
			i++
		}
	}
	accuracies = accuracies[:i]
	qualities = qualities[:i]
	rcRatios = rcRatios[:i]
	lex = lex[:i]
	trained = trained[:i]
	indices = indices[:i]

	fmt.Println("Data now over",i,"useful k-mers")

	mapSize := 50
	if k == 10 {
		mapSize = 100
	} else if k == 11 {
		mapSize = 75
	}
	rankify(accuracies,indices)
	rankify(lex,indices)
	fmt.Println("Lexicographic")
	writeHeatmap(mapSize,lex,accuracies,indices,fmt.Sprintf("%s_lex_%v.txt",alignmentsFile,k))
	rankify(qualities,indices)
	fmt.Println("quality")
	writeHeatmap(mapSize,qualities,accuracies,indices,fmt.Sprintf("%s_qual_%v.txt",alignmentsFile,k))
	rankify(rcRatios,indices)
	fmt.Println("RC balance")
	writeHeatmap(mapSize,rcRatios,accuracies,indices,fmt.Sprintf("%s_bal_%v.txt",alignmentsFile,k))
	//remove untrained from indices
	rankify(trained,indices)
	for i,v := range trained {
		if v == 0 {
			indices[i] = math.MaxInt64
		}
	}
	fmt.Println("Trained")
	writeHeatmap(mapSize,trained,accuracies,indices,fmt.Sprintf("%s_train_%v.txt",alignmentsFile,k))
	//trained, indices = rankify(trained,indices)
	//accuracy, indices = rankify(accuracy,indices)
	//writeHeatmap(50,trained,accuracies,fmt.Sprintf("%s_train_%v.txt",alignmentsFile,k))
}

//assumes getCounts has been run already
func (com *kmersCommand) rcRatios(values []float64, indices []int, k uint) {
	for i, index := range indices {
		if index == math.MaxInt64 {
			continue
		}
		rc := seeds.ReverseComplement(uint(index),k)
		forward := float64(com.kmerGoodCounts[index]+com.kmerBadCounts[index])
		backward := float64(com.kmerGoodCounts[rc]+com.kmerBadCounts[rc])
		ratio := 0.5 - forward/(forward+backward)
		if ratio < 0.0 {
			ratio = -ratio
		}
		values[i] = 0.5-ratio //range 0-0.5, with 0.5 being perfectly balanced
	}
	rankify(values, indices)
}

//assumes getCounts has been run already
func (com *kmersCommand) getQualities(values []float64, indices []int) {
	for i,index := range indices {
		if index == math.MaxInt64 {
			continue
		}
		allQ := com.kmerGoodQuality[index] + com.kmerBadQuality[index]
		//mean quality
		values[i] = float64(allQ)/float64(com.kmerGoodCounts[index]+com.kmerBadCounts[index])
	}
	rankify(values, indices)
}

func filterKmers(refSet sequence.SequenceSet, indices []int, k, numWorkers int) []int {
	//get reference kmer counts
	refKmerCounts := sequtil.KmerOccurrences(refSet.GetNSequencesFrom(0,1), k, numWorkers)
	for i, count := range refKmerCounts {
		if count == 0 {
			indices[i] = math.MaxInt64
		}
	}
	sort.Ints(indices)
	back := len(indices)-1
	for back >= 0 && indices[back] == math.MaxInt64 {
		back--
	}
	return indices[:back+1]
}

func (com *kmersCommand) doLong(numWorkers, k int, args map[string]string) {
	refSet := sequence.NewFastaSequenceSet(args["reference"], 0, 1, false, true)
	filter := false//true
	var whitelist map[int]*datum
	if filter {
		whitelist = make(map[int]*datum)
		refKmerCounts := sequtil.LongKmerOccurrences(refSet.GetNSequencesFrom(0,1), k)
		for i, count := range refKmerCounts {
			if count != 0 {
				whitelist[i] = nil
			}
		}
		refKmerCounts = nil
		runtime.GC()
		debug.FreeOSMemory()
		fmt.Println("Got whitelist from reference.")
	}

	refSeq := <-refSet.GetNSequencesFrom(0,1)
	ref := refSeq.String()
	data := getLongCounts(args["input"], args["alignment"],ref, k, whitelist)
	runtime.GC()
	debug.FreeOSMemory()
	var trainingData map[int]*datum
	if args["training"] != "" {
		tRefSet := sequence.NewFastaSequenceSet(args["training_ref"], 0, 1, false, true)
		tRefSeq := <-tRefSet.GetNSequencesFrom(0,1)
		tRef := tRefSeq.String()
		fmt.Println("Loading training data")
		trainingData = getLongCounts(args["training"], args["training_alignment"],tRef, k, data)
		fmt.Println("Loaded training data")
	} else {
		trainingData = make(map[int]*datum) //empty training data
	}
	getLongCorrelations(data, trainingData, args["alignment"], k)

	/*fmt.Println("writing kmers to",fmt.Sprintf("%s_kmers_%v.txt",args["alignment"],k))
	out, err := os.OpenFile(fmt.Sprintf("%s_kmers_%v.txt",args["alignment"],k),os.O_WRONLY | os.O_CREATE, 0755)
	if err == nil {
		for kmer, d := range data {
			if d.goodCount+d.badCount > 2 {
				out.WriteString(fmt.Sprintf("%v %v\n",sequence.KmerString(kmer,k),float64(d.goodCount)/float64(d.goodCount+d.badCount))) //just write the accuracy ranks for use as trained frequencies
			}
		}
		out.Close()
	} else {
		fmt.Println(err)
	}*/
}

func (com *kmersCommand) Run(args map[string]string) {
	numWorkers := ParseInt(args["num_workers"])
	k := ParseInt(args["k"])
	mapSize := ParseInt(args["map_size"])

	if k > 8 {
		com.doLong(numWorkers,k,args)
		return
	}
	size := 1
	for i := 0; i < k; i++ {
		size = size*4
	}
	com.kmerGoodCounts = make([]int32, size)
	com.kmerBadCounts = make([]int32, size)
	com.kmerGoodQuality = make([]int32, size)
	com.kmerBadQuality = make([]int32, size)
	indices := make([]int, size)
	for i := 0; i < len(indices); i++ {
		indices[i] = i
	}

	refSet := sequence.NewFastaSequenceSet(args["reference"], 0, 1, false, true)
	refSeq := <-refSet.GetNSequencesFrom(0,1)
	ref := refSeq.String()
	//if filtering:
	//indices = filterKmers(refSet,indices,k,numWorkers)
	refSeq = nil
	refSet = nil

	accuracies, indices := com.getCounts(args["input"], args["alignment"],ref, k, numWorkers, com.kmerGoodCounts, com.kmerBadCounts, com.kmerGoodQuality, com.kmerBadQuality, indices)

	//now indices contain only those of interest.
	//accuracies is trimmed to indices length.

	//Next stop, getting other values:
	values := make([]float64, len(indices))

	//write the lexicographic order heatmap out: accuracy vs index
	for i,index := range indices {
		values[i] = float64(index)
	}
	//values is still full-length, but now has those uninteresting values removed
	writeHeatmap(mapSize,values,accuracies,indices,fmt.Sprintf("%s_lex_%v.txt",args["alignment"],k))

	//write the quality heatmap (if quality values exist)
	if com.kmerGoodQuality != nil {
		com.getQualities(values,indices)
		writeHeatmap(mapSize,values,accuracies,indices,fmt.Sprintf("%s_qual_%v.txt",args["alignment"],k))
	}

	//write the rc-ratio heatmap
	com.rcRatios(values,indices,uint(k))
	writeHeatmap(mapSize,values, accuracies,indices,fmt.Sprintf("%s_bal_%v.txt",args["alignment"],k))

	/*
	//write out the trained heatmap (if kmer values are provided)
	if args["training"] != "" {
		_, d.others = sequtil.LoadKmerValues(args["training"])
		//loadValues(args["training"],d.others)
		for i,v := range d.others {
			if v < 2 || refKmerCounts[i] == 0 {
				d.others[i] = -2.0 //training data with no true hits is not relevant
			}
		}
		d.toRanks(true)
		writeHeatmap(mapSize,d,fmt.Sprintf("%s_train_%v.txt",args["alignment"],k))
	}
	*/
	/*out, err := os.OpenFile(fmt.Sprintf("%s_kmers_%v.txt",args["alignment"],k),os.O_WRONLY | os.O_CREATE, 0755)
	if err == nil {
		for i := 0; i < len(com.kmerGoodCounts); i++ {
			//rc := seeds.ReverseComplement(uint(i),uint(k))
			//out.WriteString(fmt.Sprintf("%v %v %v %v %v %v\n",i,com.kmerGoodCounts[i],com.kmerBadCounts[i],com.kmerGoodCounts[rc],com.kmerBadCounts[rc],d.accuracies[i]))
			//if d.others[i] > 0 && refKmerCounts[i] > 0 {//d.accuracies[i] > 0 {
			if refKmerCounts[i] > 0 && d.accuracies[i] > 0 {
				//out.WriteString(fmt.Sprintf("%v %v\n",i,d.accuracies[i])) //just write the accuracy ranks for use as trained frequencies
				out.WriteString(fmt.Sprintf("%v %v\n",sequence.KmerString(i,k),d.accuracies[i])) //just write the accuracy ranks for use as trained frequencies
				//out.WriteString(fmt.Sprintf("%v %v\n",sequence.KmerString(i,k),d.others[i]))
			}
		}
	} else {
		fmt.Println(err)
	}
	out.Close()
	*/
}

type byValue struct {
	values[]float64
	indices []int
}
func (d *byValue) Len() int {
	return len(d.indices)
}
func (d *byValue) Less(i, j int) bool {
	return d.values[i] < d.values[j]
}
func (d *byValue) Swap(i, j int) {
	d.values[i], d.values[j] = d.values[j], d.values[i]
	d.indices[i], d.indices[j] = d.indices[j], d.indices[i]
}
type byIndex struct {
	values[]float64
	indices []int
}
func (d *byIndex) Len() int {
	return len(d.indices)
}
func (d *byIndex) Less(i, j int) bool {
	return d.indices[i] < d.indices[j]
}
func (d *byIndex) Swap(i, j int) {
	d.values[i], d.values[j] = d.values[j], d.values[i]
	d.indices[i], d.indices[j] = d.indices[j], d.indices[i]
}

//rankify turns values into ranks. It assumes "ignored" values are removed and indices may not be continuous
//return max rank in values and new index length (cropping flagged indices)
func rankify(values []float64, indices []int) ([]float64,[]int) {
	d := byValue{values: values, indices:indices}
	sort.Stable(&d)
	rank := 0
	prev := -1.0
	for i, index := range indices { //now sorted by value
		if index == math.MaxInt64 {
			continue
		}
		v := values[i]
		if prev != v {
			rank++
			prev = v
		}
		values[i] = float64(rank)
	}
	fmt.Println("min-max ranks are:",values[0],rank)
	d2 := byIndex{values: values, indices:indices}
	sort.Sort(&d2)
	back := len(indices)-1
	for back >= 0 && indices[back] == math.MaxInt64 {
		back--
	}
	return values[:back+1],indices[:back+1]
}


func writeHeatmap(size int, xs []float64, ys []float64, indices []int, name string) {
	//correlation and max ranks
	maxX := 0.0
	maxY := 0.0
	meanX := 0.0
	meanY := 0.0
	count := 0
	for i, y := range ys {
		if indices[i] == math.MaxInt64 {
			continue
		}
		count++
		x := xs[i]
		meanX += x
		meanY += y
		if x > maxX {
			maxX = x
		}
		if y > maxY {
			maxY = y
		}
	}
	meanX /= float64(count)
	meanY /= float64(count)
	fmt.Println("means: ",meanX,meanY,"maxes:",maxX,maxY)
	num := 0.0
	denX := 0.0
	denY := 0.0
	for i, y := range ys {
		if indices[i] == math.MaxInt64 {
			continue
		}
		x := xs[i]
		dx := x-meanX
		dy := y-meanY
		num += dx*dy
		denX += dx*dx
		denY += dy*dy
	}
	corr := num / (math.Sqrt(denX) * math.Sqrt(denY))
	fmt.Println("Correlation: ",corr)

	hm := make([]int32, size*size)
	xRange := (float64(size-1))/maxX
	yRange := (float64(size-1))/maxY
	fmt.Println("Maxes:",maxX,maxY)
	for i, a := range ys {
		if indices[i] == math.MaxInt64 {
			continue
		}
		x := int(xs[i]*xRange + 0.5)
		y := int(a*yRange + 0.5)
		if x < size && y < size && x >= 0 && y >= 0 {
			hm[x + y*size]++
		} else{
			fmt.Println("BAD XY:",x,y,"from",xs[i],a,"using ranges",xRange,yRange)
		}
	}
	hout, err := os.OpenFile(name,os.O_WRONLY | os.O_CREATE, 0755)
	if err == nil {
		for x := 0; x < size; x++ {
			for y := 0; y < size; y++ {
				hout.WriteString(fmt.Sprintf("%v %v %v\n",x,y,hm[x+y*size]))
			}
		}
	}
	hout.Close()
}

func loadValues(filename string, values []float64) {
	infile, err := os.Open(filename)
	defer infile.Close()
	for i := 0; i < len(values); i++ {
		values[i] = -2.0
	}
	if err == nil {
		bin := bufio.NewReader(infile)
		for buf, berr := bin.ReadBytes('\n'); len(buf) > 0 || berr == nil; buf, berr = bin.ReadBytes('\n') {
			//split into two tokens
			if len(buf) == 0 {
				continue
			}
			s := string(buf[:len(buf)-1])
			tokens := strings.Split(s, " ")
			if len(tokens) > 1 {
				kmer := ParseInt(tokens[0])
				if kmer >= len(values) {
					fmt.Println("Overflow in training values.Kmer:",kmer," / ",len(values))
					continue
				}
				values[kmer] = float64(ParseInt(tokens[1]))
			}
		}
	} else {
		fmt.Println(err)
	}
}

func (com *kmersCommand) alignmentWorker(alignments <-chan *formats.SAMAlignment, seqSet sequence.SequenceSet, ids map[string]int, ref string , k int, kmerGoodCounts, kmerBadCounts, kmerGoodQuality, kmerBadQuality []int32, done chan<- bool) {
	for a := range alignments {
		seq := <-seqSet.GetNSequencesFrom(ids[a.NameA],1)
		original := seq.String()
		if a.ReverseComplement {
			seq = seq.ReverseComplement()
		}
		s := seq.String()
		q := seq.Quality()

		matches := a.AlignmentString.GetKmerMatches(k)
		prevSPos := 0
		for seqIndex, ok := <-matches; ok; seqIndex, ok = <-matches {
			refIndex, _ := <-matches
			refIndex += a.StartB

			if prevSPos == 0 { //ignore the start and end of the cigar
				prevSPos = seqIndex
			}
			sKmer := sequence.KmerValue(original[len(original)-k-seqIndex:len(original)-seqIndex]) //same position in the original
			com.lock.Lock()
			if ref[refIndex:refIndex+k] == s[seqIndex:seqIndex+k] {
				kmerGoodCounts[sKmer]++
				if q != nil && kmerGoodQuality != nil {
					kmerGoodQuality[sKmer]+=int32(q[seqIndex+k/2])
				}
			} else {
				com.kmerBadCounts[sKmer]++
				if q != nil && kmerGoodQuality != nil {
					kmerBadQuality[sKmer]+=int32(q[seqIndex+k/2])
				}
			}
			//k-mer counts:
			for prevSPos < seqIndex {
				sKmer = sequence.KmerValue(original[len(s)-k-prevSPos:len(s)-prevSPos])
				kmerBadCounts[sKmer]++
				if q != nil && kmerGoodQuality != nil {
					kmerBadQuality[sKmer]+=int32(q[prevSPos+k/2])
				}
				prevSPos++
			}
			com.lock.Unlock()
			prevSPos = seqIndex+1
		}
	}
	done <- true
}
