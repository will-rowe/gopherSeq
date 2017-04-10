/*

documentation to go here

*/
package main

import (
	"fmt"
	"io"
	"os"
	"sort"

	"github.com/will-rowe/gopherSeq/align"
	"github.com/will-rowe/gopherSeq/envtest"
	"github.com/will-rowe/gopherSeq/qcheck"
	"github.com/will-rowe/gopherSeq/version"
)

const border string = "-----------------------------------------------"

// define a struct to store information for each package in gopher-seq
type package_info struct {
	help          string
	main_function func()
}

// set up a map for all the packages in gopher-seq
var packages = map[string]package_info{
	"qcheck":  package_info{"\tquality check WGS data", qcheck.Main},
	"align":   package_info{"\talign, SNPcall and generate pseudogenome for WGS data", align.Main},
	"envtest": package_info{"\ttest runtime environment for required software", envtest.Main},
	"version": package_info{"\tprints version and exits", version.Main},
}

// create a function to print info on our packages
func printInfo() {
	var my_writer io.Writer = os.Stdout
	// Fprintf prints to an io.Writer according to format specified
	fmt.Fprintf(my_writer, "\n%s\n\t- gopherSeq -\n\na collection of tools for bacterial WGS data\n%s\n\nrun a command to find out full usage information\n\nusage:\n\tgopherSeq <command> [options]\n\ncommands:\n", border, border)
	// get all the keys from the packages map and put them in a slice of strings
	var package_names []string
	for i := range packages {
		package_names = append(package_names, i)
	}
	// sort our package names (a slice of strings) in increasing order
	sort.Strings(package_names)
	// print out package name and the help string
	for _, my_package := range package_names {
		fmt.Fprintf(my_writer, "\t%s    \t%s\n", my_package, packages[my_package].help)
	}
	fmt.Fprintf(my_writer, "\n\n")
	os.Exit(0)
}

// the main gopherSeq program gets the user command and runs the specified package
func main() {
	// print the program help if no commands were provided
	if len(os.Args) < 2 {
		printInfo()
	}
	// check that the supplied command is recognised
	var selected_package package_info
	var ok bool
	if selected_package, ok = packages[os.Args[1]]; !ok {
		fmt.Printf("unrecognised command: %s\n\n", os.Args[1])
		printInfo()
	}

	// remove the package name from the program call
	os.Args = append(os.Args[:1], os.Args[2:]...)
	selected_package.main_function()
}
