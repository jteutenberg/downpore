package main

import (
	"fmt"
	"github.com/jteutenberg/downpore/commands"
	"log"
	"os"
	"strings"
)

func alignedPrint(lines [][]string) {
	maxLengths := make([]int, 0, 10)
	for _, line := range lines {
		for i, part := range line {
			for len(maxLengths) <= i {
				maxLengths = append(maxLengths, 0)
			}
			if maxLengths[i] < len(part) {
				maxLengths[i] = len(part)
			}
		}
	}
	for _, line := range lines {
		for i, part := range line {
			fmt.Print(part)
			for j := maxLengths[i] - len(part) + 2; j > 0; j-- {
				fmt.Print(" ")
			}
		}
		fmt.Println()
	}
}

func parseArgs(com commands.Command) map[string]string {
	args, alias, _ := com.GetArgs()
	invertAlias := make(map[string]string)
	for k, v := range alias {
		invertAlias[v] = k
	}
	for i := 2; i < len(os.Args); i += 2 {
		name := strings.TrimLeft(os.Args[i], "-")
		if a, exists := invertAlias[name]; exists {
			name = a
		}
		if _, exists := args[name]; !exists {
			log.Fatal("Unrecognised argument:", name)
		}
		args[name] = os.Args[i+1]
	}
	return args
}

func main() {
	coms := []commands.Command{commands.NewOverlapCommand(), commands.NewMapCommand(), commands.NewTrimCommand(),commands.NewSubSeqCommand(), commands.NewConsensusCommand(), commands.NewAlignCommand(), commands.NewFullDenovoCommand(),commands.NewKmersCommand()}
	if len(os.Args) == 1 {
		fmt.Println("Available commands:\n help <command> Describe the command and its arguments")
		for _, com := range coms {
			fmt.Println(" " + com.GetName())
		}
	} else if os.Args[1] == "help" {
		if len(os.Args) > 2 {
			for _, com := range coms {
				if com.GetName() == os.Args[2] {
					args, alias, desc := com.GetArgs()
					lines := make([][]string, 0, len(args))
					for arg, def := range args {
						if a, exists := alias[arg]; exists {
							line := []string{"-" + arg, "-" + a, desc[arg], "(default:" + def + ")"}
							//fmt.Println("-"+arg,"\t-"+a,"\t"+desc[arg],"(default:"+def+")")
							lines = append(lines, line)
						} else {
							line := []string{"-" + arg, "", desc[arg], "(default:" + def + ")"}
							//fmt.Println("-"+arg,"\t\t"+desc[arg],"(default:"+def+")")
							lines = append(lines, line)
						}
					}
					alignedPrint(lines)
					return
				}
			}
		}
		fmt.Println("Usage: downpore help <command>\nTo see a list of available commands just run downpore")
	} else {
		for _, com := range coms {
			if com.GetName() == os.Args[1] {
				com.Run(parseArgs(com))
				return
			}
		}
		fmt.Println("Available commands:\n help <command> Describe the command and its arguments")
	}
}
