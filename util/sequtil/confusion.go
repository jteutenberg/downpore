package sequtil

import (
	"bufio"
	"github.com/jteutenberg/downpore/sequence"
	"log"
	"os"
	"strconv"
	"strings"
)

func LoadConfusionMatrix(filename string) ([][]uint8,int) {
	var matrix [][]uint8
	var k int
	fin, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer fin.Close()
	bin := bufio.NewReader(fin)
	for buf, berr := bin.ReadBytes('\n'); len(buf) > 0 || berr == nil; buf, berr = bin.ReadBytes('\n') {
		s := string(buf[:len(buf)-1])
		tokens := strings.Split(s, " ")
		if len(tokens) < 3 {
			continue
		}
		//the first sequence.. sort out storage
		if k == 0 {
			k = len(tokens[0])
			count := 1
			for i := 0; i < k; i++ {
				count *= 4
			}
			matrix = make([][]uint8, count)
			for i := 0; i < count; i++ {
				matrix[i] = make([]uint8,count)
				matrix[0][i] = 15 //maximum cost, filling first row
			}
			for i := 1; i < count; i++ {
				copy(matrix[i],matrix[0]) //fill other rows with max cost
				matrix[i][i] = 0 //perfect match
			}
			matrix[0][0] = 0
		}
		fromKmer := sequence.KmerValue(tokens[0])
		//now read values and k-mers from the remaining tokens
		for i := 1; i < len(tokens)-1; i += 2 {
			cost, _ := strconv.Atoi(tokens[i])
			kmerValue := sequence.KmerValue(tokens[i+1])
			if cost == 0 || cost > 15 {
				log.Println("Bad confusion matrix value in ",s)
				cost = 15 //an error?
			}
			matrix[fromKmer][kmerValue] = uint8(cost)
		}

	}
	return matrix,k
}
