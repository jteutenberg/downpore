package sequtil

import (
	"bufio"
	"github.com/jteutenberg/downpore/sequence"
	"log"
	"os"
	"strconv"
	"strings"
)

//LoadKmerValues reads kmer-value pairs from a file. It expects k to be
//15 or lower with values for every kmer to be present.
func LoadKmerValues(filename string) (int, []float64) {
	var values *[]float64
	var k int
	fin, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer fin.Close()

	bin := bufio.NewReader(fin)
	for buf, berr := bin.ReadBytes('\n'); len(buf) > 0 || berr == nil; buf, berr = bin.ReadBytes('\n') {
		//split into two tokens
		s := string(buf[:len(buf)-1])
		tokens := strings.Split(s, " ")
		seq := sequence.NewByteSequence(0, tokens[0], nil)
		v := seq.KmerAt(0, seq.Len())
		//the first sequence.. sort out storage
		if k == 0 {
			k = seq.Len()
			count := 1
			for i := 0; i < k; i++ {
				count *= 4
			}
			vs := make([]float64, count, count)
			values = &vs
		}
		(*values)[v], _ = strconv.ParseFloat(tokens[1], 64)
		//test for kmers that make poor seed
		if tokens[0][1:] == tokens[0][:len(tokens[0])-1] || tokens[0][2:] == tokens[0][:len(tokens[0])-2] {
			(*values)[v] = 0.0
		}
	}
	return k, *values
}
