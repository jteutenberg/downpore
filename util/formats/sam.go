package formats

import (
	"bufio"
	"os"
	"strings"
	"fmt"
)

type Cigar string

type SAMAlignment struct {
	NameA string
	NameB string
	AlignmentString Cigar
	StartA int
	StartB int
	ReverseComplement bool
}

func LoadSAM(filename string) <-chan *SAMAlignment {
	out := make(chan *SAMAlignment, 10)
	inFile, err := os.Open(filename)
	if err == nil {
		bin := bufio.NewReader(inFile)
		go func() {
			for b, err := bin.ReadBytes('\n'); err == nil || len(b) > 0; b, err = bin.ReadBytes('\n') {
				s := string(b)
				if s[0] == '@' {
					continue
				}
				tokens := strings.Fields(s)
				if len(tokens) < 6 || tokens[5] == "*" {
					continue
				}
				flags := toInt(tokens[1])
				sam := SAMAlignment{NameA:tokens[0], NameB:tokens[2], AlignmentString: Cigar(tokens[5]), StartA:0, StartB: toInt(tokens[3])-1, ReverseComplement: (flags & 0x10) != 0}
				out <- &sam
			}
			inFile.Close()
			fmt.Println(err)
			close(out)
		}()
	} else {
		fmt.Println(err)
	}
	return out
}

func (c *Cigar) CountMatches(k int) int {
	count := 0
	cigar := string(*c)
	i := 0
	for i < len(cigar) {
		j := i+1
		for cigar[j] <= '9' {
			j += 1
		}
		if cigar[j] == 'M' {
			n := toInt(cigar[i:j])
			if n >= k {
				count += n-k+1
			}
		}
		i = j+1
	}
	return count
}

//Length gets the length of sequence,reference matched
func (c *Cigar) Length() (int,int) {
	aCount := 0
	bCount := 0
	cigar := string(*c)
	i := 0
	for i < len(cigar) {
		j := i+1
		for cigar[j] <= '9' {
			j += 1
		}
		n := toInt(cigar[i:j])
		op := cigar[j]
		if op == 'M' || op == 'X' || op == '=' {
			aCount += n
			bCount += n
		} else if op == 'D' || op == 'N' {
			bCount += n
		} else if op == 'I' || op == 'H' || op == 'S' {
			aCount += n
		}
		i = j+1
	}
	return aCount,bCount
}

//GetKmerMatches retrieves pairs of int indices (query then reference) for matching
//k-mers in the cigar string
func (c *Cigar) GetKmerMatches(k int) <-chan int {
	out := make(chan int, 10)
	go func(){
		i := 0
		seq_index := 0
		ref_index := 0
		cigar := string(*c)
		for i < len(cigar) {
			j := i+1
			for cigar[j] <= '9' {
				j += 1
			}
			n := toInt(cigar[i:j])
			op := cigar[j]
			//get any matches
			if op == 'M' && n >= k {
				for m := 0; m < n-k+1; m++ {
					out <- seq_index+m
					out <- ref_index+m
				}
			}
			//move forward
			if op == 'M' || op == 'X' || op == '=' {
				seq_index += n
				ref_index += n
			} else if op == 'D' || op == 'N' {
				ref_index += n
			} else if op == 'I' || op == 'H' || op == 'S' {
				seq_index += n
			}
			i = j+1
		}
		close(out)
	}()
	return out
}
