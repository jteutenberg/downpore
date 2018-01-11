package commands

import (
	"log"
	"sort"
	"strconv"
)

type Command interface {
	GetName() string
	//GetArgs returns this commands arguments as a map of switch->default value and a map of descriptions
	GetArgs() (map[string]string, map[string]string, map[string]string)
	//Run executes the command using the given map of switch->argument value
	Run(map[string]string)
}

//MakeArgs returns an argument map filled with defaults, a mapping to aliases and a mapping to descriptions
func MakeArgs(names []string, defaults []string, descriptions []string) (map[string]string, map[string]string, map[string]string) {
	args := make(map[string]string)
	alias := make(map[string]string)
	desc := make(map[string]string)
	for i, name := range names {
		args[name] = defaults[i]
		desc[name] = descriptions[i]
	}
	//make the abreviated versions
	sort.Strings(names)
	for i := 0; i < len(names); i++ {
		if i == len(names)-1 || names[i][0] != names[i+1][0] {
			//easy case
			alias[names[i]] = string(names[i][0:1])
		} else {
			//find the minimal length to disambiguate ALL args with this first letter
			j := i + 1
			minLen := 1
			for j < len(names) && names[j][0] == names[i][0] {
				sameCount := 1
				for sameCount < len(names[j]) && sameCount < len(names[j-1]) && names[j][sameCount] == names[j-1][sameCount] {
					sameCount++
				}
				if sameCount >= minLen {
					minLen = sameCount + 1
				}
				j++
			}
			//add all those aliases
			if minLen < 4 {
				for _, n := range names[i:j] {
					alias[n] = string(n[:minLen])
				}
			}
			i = j - 1
		}
	}
	return args, alias, desc
}

func ParseInt(arg string) int {
	if s, err := strconv.ParseInt(arg, 10, 32); err == nil {
		return int(s)
	}
	log.Fatal("Invalid integer argument value:", arg)
	return 0
}
func ParseFloat(arg string) float64 {
	if s, err := strconv.ParseFloat(arg, 64); err == nil {
		return s
	}
	log.Fatal("Invalid float argument value:", arg)
	return 0
}
func ParseBool(arg string) bool {
	return arg == "1" || (len(arg) > 0 && (arg[0] == 'T' || arg[0] == 'T'))
}
