package main

import (
	"fmt"
	"github.com/jteutenberg/downpore/commands"
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
			fmt.Fprint(os.Stderr, part)
			for j := maxLengths[i] - len(part) + 2; j > 0; j-- {
				fmt.Fprint(os.Stderr, " ")
			}
		}
		fmt.Fprintln(os.Stderr, "")
	}
}

func parseArgs(com commands.Command) (map[string]string, bool) {
	args, alias, _ := com.GetArgs()
	invertAlias := make(map[string]string)
	for k, v := range alias {
		invertAlias[v] = k
	}
	for i := 2; i < len(os.Args); i++ {
		name := strings.TrimLeft(os.Args[i], "-")
		//check for = style
		equals := strings.Index(name, "=")
		var value string
		if equals > 0 {
			value = name[equals+1:]
			name = name[:equals]
		}
		if a, exists := invertAlias[name]; exists {
			name = a
		}
		if _, exists := args[name]; !exists {
			fmt.Fprintln(os.Stderr, "Unrecognised argument:", name)
			return nil, false
		}
		if equals > 0 {
			args[name] = value
		} else {
			if i >= len(os.Args)+1 {
				return nil, false
			}
			args[name] = os.Args[i+1]
			i++
		}
	}
	return args, true
}

func printHelp(com commands.Command) {
	args, alias, desc := com.GetArgs()
	lines := make([][]string, 0, len(args))
	for arg, def := range args {
		if a, exists := alias[arg]; exists {
			line := []string{"-" + arg, "-" + a, desc[arg], "(default:" + def + ")"}
			lines = append(lines, line)
		} else {
			line := []string{"-" + arg, "", desc[arg], "(default:" + def + ")"}
			lines = append(lines, line)
		}
	}
	alignedPrint(lines)
}

func main() {
	coms := []commands.Command{commands.NewTrimCommand()}
	if len(os.Args) == 1 {
		fmt.Fprintln(os.Stderr, "Available commands:\n help <command> Describe the command and its arguments")
		for _, com := range coms {
			fmt.Fprintln(os.Stderr, " "+com.GetName())
		}
	} else if os.Args[1] == "help" {
		if len(os.Args) > 2 {
			for _, com := range coms {
				if com.GetName() == os.Args[2] {
					printHelp(com)
					return
				}
			}
			fmt.Fprintln(os.Stderr, "Usage: downpore help <command>\nAvailable commands:\n")
			for _, com := range coms {
				fmt.Fprintln(os.Stderr, " "+com.GetName())
			}
		}
		fmt.Fprintln(os.Stderr, "Usage: downpore help <command>\nTo see a list of available commands just run downpore")
	} else {
		for _, com := range coms {
			if com.GetName() == os.Args[1] {
				if args, ok := parseArgs(com); !ok {
					printHelp(com)
				} else {
					com.Run(args)
				}
				return
			}
		}
		fmt.Fprintln(os.Stderr, "Available commands:\n help <command> Describe the command and its arguments")
		for _, com := range coms {
			fmt.Fprintln(os.Stderr, " "+com.GetName())
		}
	}
}
