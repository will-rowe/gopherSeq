/*

This package is a simple pipeline to align WGS data to a reference genome.

The steps included are:

 * collect sample information (paired/single-end etc.)
 * generate indices for a reference (bowtie2, faidx + fasta dict)
 * runs Bowtie2 alignment
 * processes alignment files
 * runs GATK indel correction
 * calls variants using mpileup and bcftools
 * creates a pseudogenome for each sample (modified reference sequence for each sample based on identified SNPs)

*/

package align


///////////////
// IMPORTS
//////////////
import (
  "time"
  "fmt"
  "os"
  "os/exec"
  "path"
  "strings"
  "io"
  "log"
  "runtime"
  "sync"
  "strconv"

  "github.com/will-rowe/gopherSeq/envtest"
  "github.com/alexflint/go-arg"
)


///////////////
// STRUCTS
//////////////
type sample_information struct {
  path_to_reads_1 string
  path_to_reads_2 string
  compressed bool
  paired bool
  path_to_bam string
}

type sample_list map[string]*sample_information


///////////////
// GLOBALS
//////////////
const border string = "-----------------------------------------------"
var stamp = time.Now().Format(time.RFC3339)
var samples = make(sample_list)
var (
  errorlog *os.File
  logger *log.Logger
  wg sync.WaitGroup
  threads string
)

// set up command line arguments
var args struct {
  Input []string `arg:"positional,help:input fastq files (can be .gz)"`
  Reference string `arg:"required,-r,help:specify a reference sequence (in fasta format)"`
  Output_dir string `arg:"-o,help:specify output directory"`
  Threads int `arg:"-t,help:number of processors to use [default: maximum]"`
  Keep bool `arg:"-k,help:keep temporary files [default: false]"`
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
  fmt.Fprintf(my_writer, "\n%s\n\t\t\t- gopherSeq -\n\na collection of tools for bacterial WGS data\n%s\n\nabout:\n\taligns sequencing data to a reference, calls SNPs & creates a pseudogenome for each sample\n\nusage:\n\tgopherSeq align [options] INPUT\n\nhelp:\n\tgopherSeq align --help\n", border, border)
  fmt.Fprintf(my_writer, "\n\n")
  os.Exit(0)
}


/*
  function to check user supplied arguments
*/
func argCheck() {
  var my_writer io.Writer = os.Stdout
  args.Output_dir = "./gopherSeq-align-" + string(stamp)

  // parse the ARGs
  arg.MustParse(&args)

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
    if ext := strings.HasSuffix(input_file, ".gz"); ext == true {
      getSampleInfo(strings.TrimSuffix(input_file, ".gz"), true)
    } else {
      getSampleInfo(input_file, false)
    }
  }

  // create the output directories
  if _, err := os.Stat(args.Output_dir); os.IsNotExist(err) {
    if err := os.Mkdir(args.Output_dir, 0700); err != nil {
      fmt.Fprintf(my_writer, "can't make output directory: %v\n", args.Output_dir)
      os.Exit(1)
    }
  }
  if err := os.Mkdir(args.Output_dir + "/tmp", 0700); err != nil {
    fmt.Fprintf(my_writer, "can't make tmp dir in output directory - already exists?\n")
    os.Exit(1)
  }
  if err := os.Mkdir(args.Output_dir + "/bams", 0700); err != nil {
    fmt.Fprintf(my_writer, "can't make bams dir in output directory - already exists?\n")
    os.Exit(1)
  }
  if err := os.Mkdir(args.Output_dir + "/bcfs", 0700); err != nil {
    fmt.Fprintf(my_writer, "can't make bcfs dir in output directory - already exists?\n")
    os.Exit(1)
  }
  if err := os.Mkdir(args.Output_dir + "/pseudogenomes", 0700); err != nil {
    fmt.Fprintf(my_writer, "can't make pseudogenomes dir in output directory - already exists?\n")
    os.Exit(1)
  }

  // set number of threads to use
  if args.Threads <= 0 || args.Threads > runtime.NumCPU() {
    threads = strconv.Itoa(runtime.NumCPU())
  } else {
    threads = strconv.Itoa(args.Threads)
  }
}


/*
  function to collect information from each input fastq file
*/
func getSampleInfo(input_file string, compressed bool) {
  var sample string
  var paired bool
  var path_to_reads string

  // check for .fq/.fastq extension and remove it
  if ext := strings.HasSuffix(input_file, ".fastq"); ext == true {
    sample = strings.TrimSuffix(path.Base(input_file), ".fastq")
  } else if ext := strings.HasSuffix(input_file, ".fq"); ext == true {
    sample = strings.TrimSuffix(path.Base(input_file), ".fq")
  } else {
    fmt.Printf("file format not supported: %v\n", input_file)
    os.Exit(1)
  }

  // gather sample info
  if ext := strings.HasSuffix(sample, "_1"); ext == true {
    paired = true
    if compressed == true {
      path_to_reads = input_file + ".gz" } else { path_to_reads = input_file
    }
    sample = strings.TrimSuffix(sample, "_1")
  } else if ext := strings.HasSuffix(sample, "_2"); ext == true {
    paired = true
    if compressed == true {
      path_to_reads = input_file + ".gz" } else { path_to_reads = input_file
    }
    sample = strings.TrimSuffix(sample, "_2")
  } else {
    paired = false
    if compressed == true {
      path_to_reads = input_file + ".gz" } else { path_to_reads = input_file
    }
  }

  // save sample information or append if sample basename already exists
  if _, ok := samples[sample]; ok != true {
    samples[sample] = &sample_information{path_to_reads, "", compressed, paired, ""}
  } else {
    samples[sample].path_to_reads_2 = path_to_reads
  }
}


/*
  function to set up logging
*/
func getLogging() {
  errorlog, err := os.OpenFile(args.Output_dir + "/log.txt",  os.O_RDWR|os.O_CREATE|os.O_APPEND, 0700)
  if err != nil {
    fmt.Printf("error opening log file: %v", err)
    os.Exit(1)
  }
  logger = log.New(errorlog, "gopherSeq: ", log.Lshortfile|log.LstdFlags)
  logger.Printf("--- started gopherSeq align ---\n")
}


/*
  function to generate indices from reference
*/
func createIndex() {
  bt2_reference := args.Output_dir + "/tmp/bt2_ref"
  if err := exec.Command("bowtie2-build", args.Reference, bt2_reference).Run(); err != nil {
    logger.Printf(" * couldn't create the bowtie2 index! Check the reference sequence\n")
    os.Exit(1)
  }
  index_cmd := "cp " + args.Reference + " " + args.Output_dir + "/tmp/reference.fa && samtools faidx " + args.Output_dir + "/tmp/reference.fa && samtools dict " + args.Output_dir + "/tmp/reference.fa >" + args.Output_dir + "/tmp/reference.dict"
  if err := exec.Command("bash", "-c", index_cmd).Run(); err != nil {
    fmt.Println(index_cmd)
    logger.Printf(" * couldn't create the faidx index! Check the reference sequence\n")
    os.Exit(1)
  }
}


/*
  function to run Bowtie2
*/
func runBowtie2() {
  bt2_reference := args.Output_dir + "/tmp/bt2_ref"

  // loop through samples, running one alignment at a time
  for sample, info := range samples {
    logger.Printf(" * aligning reads from %s", sample)
    outfile := args.Output_dir + "/tmp/alignment_file." + sample + ".sorted.bam"
    bt2_cmd := "bowtie2 -x " + bt2_reference + " --rg-id foo --rg SM:bar’ -p " + threads

    // customise bowtie2 command based on sample type
    if info.paired == true {
      bt2_cmd = bt2_cmd + " -1 " + info.path_to_reads_1 + " -2 " + info.path_to_reads_2
    } else {
      bt2_cmd = bt2_cmd + " -U " + info.path_to_reads_1
    }
    piped_cmd := bt2_cmd + " | samtools view -@ " + threads + " -q 10 -bh - | samtools sort -@ " + threads + " - -o " + outfile
    if err := exec.Command("bash", "-c", piped_cmd).Run(); err != nil {
      logger.Printf("failed to execute alignment: %s", piped_cmd)
      os.Exit(1)
    }

    // update sample info with bam file
    samples[sample].path_to_bam = outfile
  }
}


/*
  function to perform InDel realignment
*/
func runGATK(sample string, info *sample_information, worker int) {
  // remove duplication and index bam
  logger.Printf("\t[ worker %d: * removing duplicates and indexing %s ]", worker, sample)
  bam_nodup := args.Output_dir + "/tmp/alignment_file." + sample + ".sorted.nodup.bam"
  RMDUP := "samtools rmdup -S " + info.path_to_bam + " " + bam_nodup +  " && samtools index " + bam_nodup
  if err := exec.Command("bash", "-c", RMDUP).Run(); err != nil {
    logger.Printf("failed to execute rmdup: %s", RMDUP)
    os.Exit(1)
  }

  // create targets
  logger.Printf("\t[ worker %d: * creating targets for %s ]", worker, bam_nodup)
  RTC := "java -jar $gopherSeq_bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R " + args.Output_dir + "/tmp/reference.fa -I " + bam_nodup + " -o " + args.Output_dir + "/tmp/realigner.intervals"
  if err := exec.Command("bash", "-c", RTC).Run(); err != nil {
    logger.Printf("failed to execute create targets: %s", RTC)
    os.Exit(1)
  }

  // realign indels
  logger.Printf("\t[ worker %d: * realigning indels for %s ]", worker, bam_nodup)
  outfile := args.Output_dir + "/bams/alignment_file." + sample + ".sorted.nodup.indels_corrected.bam"
  IR := "java -jar $gopherSeq_bin/GenomeAnalysisTK.jar -T IndelRealigner -R " + args.Output_dir + "/tmp/reference.fa -I " + bam_nodup + " -targetIntervals " + args.Output_dir + "/tmp/realigner.intervals -o " + outfile
  if err := exec.Command("bash", "-c", IR).Run(); err != nil {
    logger.Printf("failed to execute indel realignment: %s", IR)
    os.Exit(1)
  }

    // update sample info with corrected bam file
    samples[sample].path_to_bam = outfile
}


/*
  function to run variant call against reference
*/
func runSNPcall(sample string, info *sample_information, worker int) {
  // run mpileup
  logger.Printf("\t[ worker %d: * running mpileup on %s ]", worker, sample)
  MPILEUP := "samtools mpileup -d 1000 -DSugBf " + args.Output_dir + "/tmp/reference.fa " + info.path_to_bam + " > " + args.Output_dir + "/tmp/" + sample + ".tmp.bcf"
  if err := exec.Command("bash", "-c", MPILEUP).Run(); err != nil {
    logger.Printf("failed to run mpileup: %s", MPILEUP)
    os.Exit(1)
  }

  // run bcftools
  logger.Printf("\t[ worker %d: * running bcftools on %s ]", worker, sample)
  outfile := args.Output_dir + "/bcfs/" + sample + ".bcf"
  BCFTOOLS := "bcftools call -c --ploidy 1 " + args.Output_dir + "/tmp/" + sample + ".tmp.bcf -O u -o " + outfile
  if err := exec.Command("bash", "-c", BCFTOOLS).Run(); err != nil {
    logger.Printf("failed to run bcftools:%s", BCFTOOLS)
    os.Exit(1)
  }

  // create pseudogenome
  logger.Printf("\t[ worker %d: * creating pseudogenome for %s ]", worker, sample)
  pseudogenome := args.Output_dir + "/pseudogenomes/" + sample + ".pseudogenome.fa"
  PSEUDO := "bcftools view " + outfile + " | $gopherSeq_bin/vcfutils.pl vcf2fa -d 5 > " + pseudogenome
  if err := exec.Command("bash", "-c", PSEUDO).Run(); err != nil {
    logger.Printf("failed to generate pseudogenome: %s", PSEUDO)
    os.Exit(1)
  }
}


/*
  function to submit tasks to worker goroutines
*/
func runGoroutines() {
  logger.Printf(" * launching goroutines")

  // setup the variables
  numberGoroutines, _ := strconv.Atoi(threads)
  taskLoad := len(samples)

  // create a buffered channel to manage the task load
  tasks := make(chan string, taskLoad)

  // launch goroutines to handle tasks
  wg.Add(numberGoroutines)
  for gr := 1; gr <= numberGoroutines; gr++ {
    go worker(tasks, gr)
  }

  // add the work to the task list
  for sample := range samples {
    tasks <- sample
  }

  // close the channel so all the Goroutines will terminate when all tasks are done
  close(tasks)

  // wait for the tasks to be completed
  wg.Wait()
}


/*
  function to run worker goroutines to complete tasks from a list
*/
func worker(tasks chan string, worker int) {
  // send completion signal
  defer wg.Done()

  // endless for loop to receive tasks (until channel closed and empty)
  for {
    // get a sample from the task list
    sample, ok := <- tasks

    // check if channel is closed
    if !ok {
      logger.Printf("\t[ worker %d: shutting down ]", worker)
      return
    }

    // otherwise, start work
    logger.Printf("\t[ worker %d: starting task ]", worker)
    info := samples[sample]

    // run the task
    runGATK(sample, info, worker) // indel correction
    runSNPcall(sample, info, worker)  // variant call

    // notify task completion
    logger.Printf("\t[ worker %d: completed task ]", worker)
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
    argCheck()
  }

  // spinner (this was just to play with goroutines)
  go Spinner(100*time.Millisecond)

  // start the logger
  getLogging()
  defer errorlog.Close()

  // check for gopherSeq bin
  logger.Printf("checking for gopherSeq bin . . .")
  passed, messages := envtest.BinCheck()
  for _, message := range messages {
    logger.Printf("%v", message)
  }
  if passed == false {
    logger.Printf("gopherSeq bin check failed!\n")
    os.Exit(1)
  }

  // check for required programs
  logger.Printf("checking for required software . . .")
  passed, messages = true, nil
  passed, messages = envtest.Test4align_progs()
  for _, message := range messages {
    logger.Printf("%v", message)
  }
  if passed == false {
    logger.Printf("program check failed!\n")
    os.Exit(1)
  }

  // print some messages
  logger.Printf("checking for input arguments . . .")
  logger.Printf(" * reference sequence supplied --> %v", args.Reference)
  logger.Printf(" * number of samples supplied --> %d", len(samples))
  for sample, information := range samples {
    logger.Printf("\tSAMPLE=%v READS=%v %v", sample, information.path_to_reads_1, information.path_to_reads_2)
    // include a check for paired end samples
    if information.paired == true {
      if len(information.path_to_reads_1) == 0 || len(information.path_to_reads_2) == 0 {
        logger.Printf("\tonly one read file found for this sample - the pipeline thinks it should be paired")
        os.Exit(1)
      }
    }
  }
  logger.Printf(" * number of threads to be used --> %s", threads)
  logger.Printf(" * keeping temporary files --> %t", args.Keep)
  logger.Printf(" * output directory --> %s", args.Output_dir)

  // create Bowtie2 index
  logger.Printf("building Bowtie2 index . . .")
  createIndex()

  // run Bowtie2
  logger.Printf("--- started read alignment ---")
  logger.Printf("running Bowtie2 and sorting with Samtools . . .")
  runBowtie2()

  // run InDel correction and perform mpileup
  logger.Printf("--- started InDel correction & SNP call ---")
  logger.Printf("running GATK, samtools + bcftools . . .")
  runGoroutines()

  // clean up
  logger.Printf("--- finished ---")
  logger.Printf("all files saved to: %s", args.Output_dir)
  if args.Keep == false {
    if err := os.RemoveAll(args.Output_dir + "/tmp"); err != nil {
      logger.Printf("could not remove the tmp file directory!")
    } else {
      logger.Printf("removed temporary files")
    }
  }
}
