/*

This package is a simple QC pipeline to check WGS data

The steps included are:

 *
 * launch the align pipeline (optional)

*/

package qcheck

///////////////
// IMPORTS
//////////////
import (
	"fmt"
	"io"
	"os"
	"os/exec"
	"path"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/alexflint/go-arg"
	"github.com/will-rowe/gopherSeq/envtest"
)

///////////////
// GLOBALS
//////////////
const border string = "-----------------------------------------------"

var stamp = time.Now().Format(time.RFC3339)
var threads string
var trimmed_samples []string

// set up command line arguments
var args struct {
	Input      []string `arg:"positional"`
	Output_dir string   `arg:"-o,help:specify output directory "`
	Threads    int      `arg:"-t,help:number of processors to use [default: maximum]"`
	Align      bool     `arg:"-a,help:run align pipeline after the QC check finishes [default: false]"`
	Reference  string   `arg:"-r,help:specify a reference sequence (required if --align selected)"`
}

///////////////
// FUNCTIONS
//////////////
/*
  function to print info on our package
*/
func printInfo() {
	var my_writer io.Writer = os.Stdout

	// Fprintf prints to an io.Writer according to format specified
	fmt.Fprintf(my_writer, "\n%s\n\t- gopherSeq -\n\na collection of tools for bacterial WGS data\n%s\n\nabout:\n\truns basic QC on fastq files\n\nusage:\n\tgopherSeq qcheck [options] INPUT\n\nhelp:\n\tgopherSeq qcheck --help\n", border, border)
	fmt.Fprintf(my_writer, "\n\n")
	os.Exit(0)
}

/*
  function to check user supplied arguments
*/
func argCheck() {
	var my_writer io.Writer = os.Stdout
	args.Output_dir = "./gopherSeq-qcheck-" + string(stamp)

	// parse the ARGs
	arg.MustParse(&args)

	// check the input files exist
	for _, input_file := range args.Input {
		if _, err := os.Stat(input_file); err != nil {
			if os.IsNotExist(err) {
				fmt.Fprintf(my_writer, "file does not exist: %v\n", input_file)
				os.Exit(1)
			} else {
				fmt.Fprintf(my_writer, "can't access file: %v\n", input_file)
				os.Exit(1)
			}
		}

		// check the file extension and then parse the filenames
		if strings.HasSuffix(input_file, ".gz") {
			input_file = strings.TrimSuffix(input_file, ".gz")
		}
		if strings.HasSuffix(input_file, ".fq") || strings.HasSuffix(input_file, ".fastq") {
			fmt.Printf(" * found input file --> %v\n", input_file)
		} else {
			fmt.Printf("a supplied input file does not seem to be in fastq format: %v\n", input_file)
			os.Exit(1)
		}
	}

	// set number of threads to use
	if args.Threads <= 0 || args.Threads > runtime.NumCPU() {
		threads = strconv.Itoa(runtime.NumCPU())
	} else {
		threads = strconv.Itoa(args.Threads)
	}

	// if align selected, make sure reference supplied
	if args.Align == true {
		// check the reference sequence exists
		if _, err := os.Stat(args.Reference); err != nil {
			if os.IsNotExist(err) {
				if os.IsNotExist(err) {
					fmt.Fprintf(my_writer, "file does not exist: %v\n", args.Reference)
					os.Exit(1)
				} else {
					fmt.Fprintf(my_writer, "can't access file: %v\n", args.Reference)
					os.Exit(1)
				}
			}
		}
	}

	// create the output directories
	if _, err := os.Stat(args.Output_dir); os.IsNotExist(err) {
		if err := os.Mkdir(args.Output_dir, 0700); err != nil {
			fmt.Fprintf(my_writer, "can't make output directory: %v\n", args.Output_dir)
			os.Exit(1)
		}
	}
	if err := os.Mkdir(args.Output_dir+"/QC_files", 0700); err != nil {
		fmt.Fprintf(my_writer, "can't make QC directory in %v\n", args.Output_dir)
		os.Exit(1)
	}

}

/*
  function to run QC
*/
func qcData() {
	gopherSeq_bin := os.Getenv("gopherSeq_bin")

	// loop through samples and run each qc program
	for _, sample := range args.Input {
		fmt.Printf("[ current sample: %v ]\n", sample)
		basename := path.Base(sample)

		// fastqc
		fmt.Println(" * running fastqc")
		fastqc_cmd := "fastqc --threads " + threads + " --quiet --outdir " + args.Output_dir + "/QC_files " + sample
		if err := exec.Command("bash", "-c", fastqc_cmd).Run(); err != nil {
			fmt.Printf("fastqc command failed: %v\n", fastqc_cmd)
			os.Exit(1)
		}

		// kraken
		fmt.Println(" * running kraken")
		if _, err := os.Stat(gopherSeq_bin + "/kraken_db"); os.IsNotExist(err) {
			fmt.Println("\t- can't find kraken_db (needs symoblic link in the gopherSeq_bin)")
			fmt.Println("\t- skipping kraken")
		} else {
			kraken_cmd := "kraken --threads " + threads + " --preload --fastq-input --gzip-compressed --db $gopherSeq_bin/kraken_db " + sample + " | kraken-report --db $gopherSeq_bin/kraken_db > " + args.Output_dir + "/QC_files/krakenreport.txt"
			if err := exec.Command("bash", "-c", kraken_cmd).Run(); err != nil {
				fmt.Printf("kraken command failed: %v\n", kraken_cmd)
				os.Exit(1)
			}
		}

		// trimmomatic
		fmt.Println(" * running trimmomatic")
		var trim_cmd string
		if _, err := os.Stat(gopherSeq_bin + "/adapters.fa"); os.IsNotExist(err) {
			fmt.Println("\t- no adapter file supplied (needs symoblic link in the gopherSeq_bin)")
			fmt.Println("\t- just performing quailty-based trimming")
			trim_cmd = "trimmomatic SE -threads " + threads + " " + sample + " " + args.Output_dir + "/QC_files/trimmed." + basename + " SLIDINGWINDOW:4:20 MINLEN:100 &> " + args.Output_dir + "/QC_files/trimmomatic_logfile_for_" + basename + ".log"
		} else {
			trim_cmd = "trimmomatic SE -threads " + threads + " " + sample + " " + args.Output_dir + "/QC_files/trimmed." + basename + " ILLUMINACLIP:$gopherSeq_bin/adapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:100 &> " + args.Output_dir + "/QC_files/trimmomatic_logfile_for_" + basename + ".log"
		}
		if err := exec.Command("bash", "-c", trim_cmd).Run(); err != nil {
			fmt.Printf("trimmomatic command failed: %v\n", trim_cmd)
			os.Exit(1)
		}

		// add timmed sample to an array (to submit to align program if asked)
		trimmed_samples = append(trimmed_samples, args.Output_dir+"/QC_files/trimmed."+basename)
	}

	// run multiqc once all samples have been run through the programs
	fmt.Println(" * running multiqc")
	multiqc_cmd := "multiqc -o " + args.Output_dir + " " + args.Output_dir
	if err := exec.Command("bash", "-c", multiqc_cmd).Run(); err != nil {
		fmt.Printf("multiqc command failed: %v\n", multiqc_cmd)
		fmt.Printf("will continue with pipeline but no multiqc report will be available\nrecommend you check your multiqc install...\n")
	}
}

/*
  function for silly spinner
*/
func Spinner(delay time.Duration) {
	for {
		for _, r := range `-\|/` {
			fmt.Printf("\r%c", r)
			time.Sleep(delay)
		}
	}
}

///////////////
// MAIN
//////////////
func Main() {
	// print usage or parse arguments
	if len(os.Args) < 2 {
		printInfo()
	} else {
		fmt.Println("starting QC check . . .")
		argCheck()
	}

	// check for gopherSeq bin
	fmt.Printf("checking for gopherSeq bin . . .\n")
	passed, messages := envtest.BinCheck()
	for _, message := range messages {
		fmt.Printf("%v", message)
	}
	if passed == false {
		fmt.Printf("\ngopherSeq bin check failed!\n")
		os.Exit(1)
	}

	// check for required programs
	fmt.Printf("checking for required software . . .\n")
	passed, messages = true, nil
	passed, messages = envtest.Test4qcheck_progs()
	for _, message := range messages {
		fmt.Printf("%v", message)
	}
	if passed == false {
		fmt.Printf("\nprogram check failed!\n")
		os.Exit(1)
	}

	// perform QC
	fmt.Println("running QC programs . . .")
	qcData()
	fmt.Println("QC finished!")

	// run the align pipeline if requested
	if args.Align == true {
		fmt.Println("now starting align pipeline on the trimmed samples . . .")

		// set up the align command
		align := []string{"align", "-t", threads, "-r", args.Reference, "-o", args.Output_dir}
		align = append(align, trimmed_samples...)

		// launch command using a goroutine and create a waitgroup to wait for completion
		var wg sync.WaitGroup
		wg.Add(1)
		go func() {
			defer wg.Done()
			if err := exec.Command("gopherSeq", align...).Run(); err != nil {
				fmt.Printf("could not run align pipeline: gopherSeq %s\n", align)
			}
		}()

		// put a spinner on the screen so user knows other pipeline is running
		go Spinner(100 * time.Millisecond)

		// wait for goroutine
		wg.Wait()
	}

}
