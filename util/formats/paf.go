package formats

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

//Overlap is one entry of a pairwise alignment file. This could be generated
//by BAM/SAM, but for now we only have a PAF reader (i.e. for minimap2)
type Overlap struct {
	NameA string //Query
	NameB string //Target
	LengthA int
	LengthB int
	StartA int
	StartB int
	EndA int
	EndB int
	ReverseComplement bool
	Matches int
	Length int
	Quality int
}

func toInt(s string) int {
	v, _ := strconv.ParseInt(s,10,64)
	return int(v) //errors just set to zero
}

func LoadPAF(filename string) <-chan *Overlap {
	out := make(chan *Overlap, 10)
	inFile, err := os.Open(filename)
	if err == nil {
		bin := bufio.NewReader(inFile)
		go func() {
			for b, err := bin.ReadBytes('\n'); err == nil; b, err = bin.ReadBytes('\n') {
				tokens := strings.Fields(string(b))
				if len(tokens) < 12 {
					continue
				}
				over := Overlap{NameA:tokens[0], NameB:tokens[5], LengthA: toInt(tokens[1]), LengthB:toInt(tokens[6]), StartA:toInt(tokens[2]), EndA: toInt(tokens[3]), StartB: toInt(tokens[7]), EndB:toInt(tokens[8]), ReverseComplement: tokens[4] == "-", Matches:toInt(tokens[9]), Length:toInt(tokens[10]), Quality:toInt(tokens[11])}
				out <- &over
			}
			inFile.Close()
			close(out)
		}()
	}
	return out
}

func (v *Overlap) StringPAF() string {
	return fmt.Sprintf("")
}

func (v *Overlap) StringGFA() string {
	return fmt.Sprintf("")
}
