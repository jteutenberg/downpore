package commands

import (
	"fmt"
)

type versionCommand string

func NewVersionCommand() Command {
	v := versionCommand("0.2")
	return &v
}

func (com *versionCommand) GetName() string {
	return "version"
}
func (com *versionCommand) GetArgs() (map[string]string, map[string]string, map[string]string) {
	s := make(map[string]string)
	return s, s, s
}
func (com *versionCommand) Run(args map[string]string) {
	fmt.Println("downpore version", *com)
}
