package sequence

import (
	"bufio"
	"fmt"
	"github.com/jteutenberg/downpore/util"
	"io"
	"log"
	"math"
	"os"
	"strings"
	"sync"
)

type ReadSeekCloser interface {
	io.ReadSeeker
	io.Closer
}

type SequenceSet interface {
	GetSequences() <-chan Sequence
	GetNSequencesFrom(int, int) <-chan Sequence //starting from previously seen sequence id
	GetSequencesByID([]int) <-chan Sequence
	GetLength(int) int //gets the length of a (previously read) sequence
	GetBases() int64   //all bases in this set (read so far)
	GetName(int) string
	SetName(int, string)
	SetIgnore(int, bool) //a sequence id that will be skipped on future requests
	SetFrontTrim(int, int)
	SetBackTrim(int, int)
	GetFrontTrim(int) int
	GetBackTrim(int) int
	Size() int
	AddSequence(Sequence, string) //additional sequences to be kept in memory, these appear after the fasta on read
	Write(io.Writer, bool)
}

type fastaSequenceSet struct {
	filename   string
	offsets    []int64 //byte offset to start of each sequence, including the frontTrim
	ignore     []bool
	frontTrim  []int
	backTrim   []int
	lengths    []int
	bases      int64
	names      []string
	cached     []Sequence
	minLen     int
	isFastq    bool //whether quality scores are available or not
	size       int
	extras     []Sequence
	extraNames []string
	numWorkers int
	cache      bool //whether or not to cache sequences in memory
	cacheFull  bool //whether the whole input file has been cached
}

func NewFastaSequenceSet(filename string, minLength int, numWorkers int, cache bool) SequenceSet {
	f := fastaSequenceSet{filename: filename, offsets: make([]int64, 0, 500000), ignore: make([]bool, 0, 500000), frontTrim: make([]int, 0, 500000), backTrim: make([]int, 0, 500000), lengths: make([]int, 0, 500000), bases: 0, names: make([]string, 0, 500000), minLen: minLength, isFastq: false, extras: make([]Sequence, 0, 20), extraNames: make([]string, 0, 20), numWorkers: numWorkers, cache: cache}
	if cache {
		f.cached = make([]Sequence, 0, 500000)
	}
	return &f
}

func (f *fastaSequenceSet) sendExtras(sentCount, maxSeqs, nextID int, seqOut chan Sequence) bool {
	if sentCount >= maxSeqs {
		return false
	}
	for _, seq := range f.extras {
		seq.setID(nextID)
		if len(f.ignore) <= nextID {
			f.ignore = append(f.ignore, false)
		}
		if !f.ignore[nextID] {
			seqOut <- seq
			sentCount++
		}
		nextID++
		if sentCount >= maxSeqs {
			close(seqOut)
			return true
		}
	}
	return false
}

func (f *fastaSequenceSet) readFasta(in ReadSeekCloser, nextID int, maxSeqs int) <-chan Sequence {
	seqOut := make(chan Sequence, 10)
	go func() {
		sentCount := 0
		var offset int64
		if f.cache && len(f.cached) > 0 {
			//send any cached through first, apply trim
			for ; nextID < len(f.cached) && sentCount < maxSeqs; nextID++ {
				if !f.ignore[nextID] {
					seqOut <- f.cached[nextID].SubSequence(f.frontTrim[nextID], f.cached[nextID].Len()-f.backTrim[nextID])
					sentCount++
				}
			}
			if !f.cacheFull {
				//need to seek ahead
				if nextID < len(f.offsets) {
					offset = f.offsets[nextID]
				} else if sentCount > 0 {
					offset = f.offsets[len(f.offsets)-1] + int64(f.cached[len(f.cached)-1].Len()+1)
				}
				if _, err := in.Seek(offset, io.SeekStart); err != nil {
					log.Fatal(err, "during seek. Stopping sequence input here.")
				}
			}

		}
		a := byte('A')
		t := byte('T')
		fastqComment := byte('@')
		plus := byte('+')
		//read any we've seen before, seeking through the file
		buf := make([]byte, 1000000, 1000000) //up to one mega-base
		for ; nextID < len(f.ignore) && sentCount < maxSeqs; nextID++ {
			if f.ignore[nextID] {
				if nextID == len(f.ignore)-1 {
					offset = f.offsets[nextID]
					if _, err := in.Seek(offset, io.SeekStart); err != nil {
						log.Println(err, "during seek. Stopping sequence input here.")
						break
					}
				}
				continue
			}
			offset = f.offsets[nextID]
			if _, err := in.Seek(offset, io.SeekStart); err != nil {
				log.Println(err, "during seek. Stopping sequence input here.")
				break
			}
			if len(buf) <= f.lengths[nextID] {
				buf = make([]byte, f.lengths[nextID]+10000, f.lengths[nextID]+10000)
			}
			if n, _ := io.ReadFull(in, buf[:f.lengths[nextID]]); n == f.lengths[nextID] {
				//transform into 0-4 values. A=65, C=67, G=71, T=84
				data := make([]byte, n, n)
				for i, b := range buf[:n] {
					data[i] = ((b >> 1) ^ ((b & 4) >> 2)) & 3
				}
				seq := byteSequence{data: data, id: nextID, name: &(f.names[nextID])}
				offset += int64(n)
				if f.isFastq { //get corresponding quality data
					quality := make([]byte, n, n)
					//skip ahead by back trim + 3 (2 endlines and a "+") + front trim
					offset += int64(f.backTrim[nextID] + f.frontTrim[nextID] + 3)
					if _, err := in.Seek(offset, io.SeekStart); err != nil {
						log.Println(err, "during seek to quality values. Stopping sequence input here.")
						break
					}
					if m, _ := io.ReadFull(in, quality); m == n {
						for i, b := range quality {
							quality[i] = b - 33
						}
						seq.quality = quality
						offset += int64(m)
					} else {
						log.Println("Unexpected quality line size", m, "not matching expected", f.lengths[nextID], "in", f.names[nextID])
						break
					}
				}
				sentCount++
				seqOut <- &seq
			} else {
				log.Println("Unexpected read size", n, "not matching expected", f.lengths[nextID])
				break //hit the end?
			}
		}
		lastName := ""
		bin := bufio.NewReader(in)
		//read to the end of one line (either leftovers from seen sequences, or a comment line)
		if b, err := bin.ReadBytes('\n'); err != nil {
			if !f.sendExtras(sentCount, maxSeqs, nextID, seqOut) {
				close(seqOut)
			}
			in.Close()
			return
		} else if b[0] == fastqComment {
			f.isFastq = true
			lastName = string(b[1:])
			offset += int64(len(b))
		} else {
			offset += int64(len(b))
			lastName = string(b[1:])
		}

		if !f.cacheFull {
			//now read any new sequences
			for buf, err := bin.ReadBytes('\n'); sentCount < maxSeqs && (len(buf) > 0 || err == nil); buf, err = bin.ReadBytes('\n') {
				if buf[0] >= a && buf[0] <= t { // a sequence
					readSeq := len(buf) >= f.minLen
					if readSeq {
						f.ignore = append(f.ignore, false)
						f.offsets = append(f.offsets, offset)
						f.frontTrim = append(f.frontTrim, 0)
						f.backTrim = append(f.backTrim, 0)
						f.lengths = append(f.lengths, len(buf)-1)
						f.names = append(f.names, strings.Fields(lastName)[0])
						f.size++
						//transform into 0-4 values. A=65, C=67, G=71, T=84
						for i, b := range buf {
							buf[i] = ((b >> 1) ^ ((b & 4) >> 2)) & 3
						}
						seq := byteSequence{data: buf[:len(buf)-1], id: nextID, name: &(f.names[len(f.names)-1])}
						f.bases += int64(len(buf) - 1)
						nextID++
						//if fastq, read error line
						if f.isFastq {
							offset += int64(len(buf))      //the sequence
							buf, err = bin.ReadBytes('\n') //+ line
							if err != nil || buf[0] != plus {
								log.Fatal("Invalid fastq format (on + line):", string(buf))
							}
							offset += int64(len(buf))
							buf, _ = bin.ReadBytes('\n') //error line
							if len(buf) == len(seq.data)+1 {
								for i, b := range buf {
									buf[i] = b - 33
								}
								seq.quality = buf[:len(buf)-1]
							}
						}
						sentCount++
						if f.cache {
							f.cached = append(f.cached, &seq)
						}
						seqOut <- &seq
					} else if f.isFastq { //manual skip
						//manually skip error line
						offset += int64(len(buf)) //the sequence
						buf, err = bin.ReadBytes('\n')
						if err != nil || buf[0] != plus {
							log.Fatal("Invalid fastq formant (on + line):", string(buf))
						}
						offset += int64(len(buf))    //the + line
						buf, _ = bin.ReadBytes('\n') //error line
					}
				} else if buf[0] == fastqComment {
					f.isFastq = true
					lastName = string(buf[1:])
					//and skip it
				} else {
					lastName = string(buf[1:])
				}
				offset += int64(len(buf))
				if err != nil {
					break
				}
			}
		}
		f.cacheFull = f.cacheFull || (f.cache && sentCount < maxSeqs) //i.e. we finished before our limit was reached
		if !f.sendExtras(sentCount, maxSeqs, nextID, seqOut) {
			close(seqOut)
		}
		in.Close()
	}()
	return seqOut
}

func (f *fastaSequenceSet) GetNSequencesFrom(index int, n int) <-chan Sequence {
	if index != 0 && index >= len(f.offsets) {
		out := make(chan Sequence)
		close(out)
		return out
	}
	var inFile ReadSeekCloser
	inFile, err := os.Open(f.filename)
	if err != nil {
		out := make(chan Sequence)
		close(out)
		return out
	}
	//test for gzip by extension (.gz)
	if strings.HasSuffix(f.filename, ".gz") {
		inFile = util.NewSeekableGZipReader(inFile)
	}
	return f.readFasta(inFile, index, n)
}

func (f *fastaSequenceSet) GetSequences() <-chan Sequence {
	return f.GetNSequencesFrom(0, int(math.MaxInt32))
}

func (f *fastaSequenceSet) GetSequencesByID(ids []int) <-chan Sequence {
	oldIgnore := f.ignore
	f.ignore = make([]bool, len(oldIgnore), len(oldIgnore))
	for i := 0; i < len(f.ignore); i++ {
		f.ignore[i] = true
	}
	for _, id := range ids {
		f.ignore[id] = false
	}
	finalReturn := make(chan Sequence, 20)
	seqs := f.GetNSequencesFrom(0, int(math.MaxInt32)) //len(ids))
	go func() {
		sentCount := 0
		for seq := range seqs {
			sentCount++
			finalReturn <- seq
		}
		f.ignore = oldIgnore
		close(finalReturn)
	}()
	return finalReturn
}

func (f *fastaSequenceSet) GetLength(id int) int {
	return f.lengths[id]
}

func (f *fastaSequenceSet) GetBases() int64 {
	return f.bases
}

func (f *fastaSequenceSet) GetName(id int) string {
	if id >= len(f.names) {
		//test extras
		for i, seq := range f.extras {
			if seq.GetID() == id {
				return f.extraNames[i]
			}
		}
	} else {
		return f.names[id]
	}
	return fmt.Sprint(id)
}

func (f *fastaSequenceSet) SetName(id int, name string) {
	if id >= len(f.names) {
		//test extras
		for i, seq := range f.extras {
			if seq.GetID() == id {
				f.extraNames[i] = name
				return
			}
		}
	} else {
		f.names[id] = name
	}
}

func (f *fastaSequenceSet) SetIgnore(id int, ignore bool) {
	f.ignore[id] = ignore
}
func (f *fastaSequenceSet) SetFrontTrim(id, trim int) {
	f.offsets[id] += int64(trim - f.frontTrim[id])
	f.lengths[id] -= trim - f.frontTrim[id]
	f.frontTrim[id] = trim
}
func (f *fastaSequenceSet) SetBackTrim(id, trim int) {
	f.lengths[id] -= trim - f.backTrim[id]
	f.backTrim[id] = trim
}
func (f *fastaSequenceSet) GetFrontTrim(id int) int {
	return f.frontTrim[id]
}
func (f *fastaSequenceSet) GetBackTrim(id int) int {
	return f.backTrim[id]
}
func (f *fastaSequenceSet) Size() int {
	return f.size
}
func (f *fastaSequenceSet) AddSequence(seq Sequence, name string) {
	f.extras = append(f.extras, seq)
	f.extraNames = append(f.extraNames, name)
}

func (f *fastaSequenceSet) fastaWriter(seqs <-chan Sequence, out *bufio.Writer, fullNames bool, lock *sync.Mutex, done chan<- bool) {
	for s := range seqs {
		var str string
		if fullNames {
			str = fmt.Sprintf(">%s\n%s\n", f.GetName(s.GetID()), s.String())
		} else {
			str = fmt.Sprintf(">%v\n%s\n", s.GetID(), s.String())
		}
		lock.Lock()
		out.WriteString(str)
		lock.Unlock()
	}
	done <- true
}
func (f *fastaSequenceSet) fastqWriter(seqs <-chan Sequence, out *bufio.Writer, fullNames bool, lock *sync.Mutex, done chan<- bool) {
	for s := range seqs {
		var str string
		quality := s.Quality()
		for i, b := range quality {
			quality[i] = b + 33 //TODO: temporary change
		}
		if fullNames {
			str = fmt.Sprintf("@%s\n%s\n+\n%s\n", f.GetName(s.GetID()), s.String(), string(quality))
		} else {
			str = fmt.Sprintf("@%v\n%s\n+\n%s\n", s.GetID(), s.String(), string(quality))
		}
		lock.Lock()
		out.WriteString(str)
		lock.Unlock()
	}
	done <- true
}

//Write re-reads the input sequences, writing out trimmed versions of non-ignored sequences
func (f *fastaSequenceSet) Write(out io.Writer, fullNames bool) {
	seqs := f.GetSequences()
	done := make(chan bool, f.numWorkers)
	var lock sync.Mutex

	bout := bufio.NewWriter(out)
	if f.isFastq {
		for i := 0; i < f.numWorkers; i++ {
			go f.fastqWriter(seqs, bout, fullNames, &lock, done)
		}
	} else {
		for i := 0; i < f.numWorkers; i++ {
			go f.fastaWriter(seqs, bout, fullNames, &lock, done)
		}
	}

	for i := 0; i < f.numWorkers; i++ {
		<-done
	}
	bout.Flush()
}
