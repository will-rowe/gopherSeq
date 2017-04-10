/*

This package tests the environment for required software.

*/

package envtest

///////////////
// IMPORTS
//////////////
import (
	"fmt"
	"io"
	"net/http"
	"os"
	"os/exec"
	"regexp"
	"strings"

	"github.com/alexflint/go-arg"
	"github.com/mitchellh/go-homedir"
)

///////////////
// GLOBALS
//////////////
const border string = "-----------------------------------------------"

// set up command line arguments for the envtest package
var args struct {
	Run bool `arg:"-r,help:test runtime environment and exit"`
}
var check_programs = []string{
	"fastqc",
	"trimmomatic",
	"kraken",
	"kraken-report",
	"multiqc",
}
var align_programs = []string{
	"bash",
	"java",
	"bowtie2",
	"samtools",
	"bcftools",
}
var messages = []string{}
var passed bool = true

///////////////
// FUNCTIONS
//////////////
/*
  function to print info on our package
*/
func printInfo() {
	var my_writer io.Writer = os.Stdout
	fmt.Fprintf(my_writer, "\n%s\n\t- gopherSeq -\n\na collection of tools for bacterial WGS data\n%s\n\nabout:\n\ttests runtime environment for required software\n\nusage:\n\tgopherSeq envtest --run\n\nhelp:\n\tgopherSeq envtest --help\n", border, border)
	fmt.Fprintf(my_writer, "\n\n")
	os.Exit(0)
}

/*
  functions to test for gopherSeq bin
*/
func BinCheck() (passed bool, messages []string) {

	// check the env variable exists
	if len(os.Getenv("gopherSeq_bin")) == 0 {
		messages = append(messages, " * gopherSeq_bin environment variable not found!\n")
		messages = append(messages, " * setting one up (~/.gopherSeq_bin)\n")

		// if it doesn't, download the files and set up the env variable
		if passed = NewBin(); passed != true {
			messages = append(messages, " * failed to set up the gopherSeq_bin!\n")
			return

			// now we have a bin, start the check again
		} else {
			messages = append(messages, " * downloaded the files and created the gopherSeq_bin!\n")
			messages = append(messages, " * please refresh your session with `source ~/.profile` and then re-try gopheSeq\n")
			passed = false
			return
		}

		// if env variable exists, test the bin for GATK, picard and vcfutils
	} else {
		passed = true
		messages = append(messages, " * found gopherSeq_bin --> "+os.Getenv("gopherSeq_bin")+"\n")
		gatk := "java -jar $gopherSeq_bin/GenomeAnalysisTK.jar -h"
		err := exec.Command("bash", "-c", gatk).Run()
		if err != nil {
			messages = append(messages, " * GATK not working - check the java and GATK install\n")
			passed = false
		} else {
			messages = append(messages, " * found GATK --> "+os.Getenv("gopherSeq_bin")+"\n")
		}
		picard := "java -jar $gopherSeq_bin/picard.jar MarkDuplicates --version"
		output, err := exec.Command("bash", "-c", picard).CombinedOutput()
		if match, _ := regexp.MatchString("2.9.0", string(output)); match == false {
			messages = append(messages, " * Picard not working - check the java and picard install\n")
			passed = false
		} else {
			messages = append(messages, " * found Picard --> "+os.Getenv("gopherSeq_bin")+"\n")
		}

		// make sure they are using the vcfutils with vcf2fa
		vcfutils := "grep vcf2fa $(which $gopherSeq_bin/vcfutils.pl)"
		response, err := exec.Command("bash", "-c", vcfutils).Output()
		if err != nil {
			messages = append(messages, " * can't find vcfutils.pl in gopherSeq bin")
			passed = false
		} else if match, _ := regexp.MatchString("vcfutils.pl vcf2fa", string(response)); match == false {
			messages = append(messages, " * your version of vcfutils.pl seems incorrect - please check for vcf2fa\n")
			passed = false
		} else {
			messages = append(messages, " * found vcfutils.pl --> "+os.Getenv("gopherSeq_bin")+"\n")
		}
	}
	return
}

/*
  function to create gopherSeq bin
*/
func NewBin() (passed bool) {

	// make new bin
	homeDir, _ := homedir.Dir()
	newBin := homeDir + "/.gopherSeq_bin/"
	if err := os.Mkdir(newBin, 0777); err != nil {
		fmt.Println("can't make gopherSeq_bin - does it already exist?\n")
		fmt.Println(err)
		passed = false
	}

	// get set of urls to download
	var urls = []string{
		"https://github.com/will-rowe/gopherSeq/raw/master/bin/GenomeAnalysisTK.jar",
		"https://github.com/will-rowe/gopherSeq/raw/master/bin/picard.jar",
		"https://github.com/will-rowe/gopherSeq/raw/master/bin/vcfutils.pl",
	}

	// download to bin
	for _, url := range urls {
		tokens := strings.Split(url, "/")
		fileName := newBin + tokens[len(tokens)-1]
		//messages = append(messages, "Downloading " + url + " to " + fileName + "\n")
		output, err := os.Create(fileName)
		if err != nil {
			//messages = append(messages, "Error while creating" + fileName + "\n")
			passed = false
			break
		}
		defer output.Close()
		response, err := http.Get(url)
		if err != nil {
			//messages = append(messages, "Error while downloading" + url + "\n")
			passed = false
			break
		}
		defer response.Body.Close()
		_, err = io.Copy(output, response.Body)
		if err != nil {
			//messages = append(messages, "Error while downloading" + url + "\n")
			passed = false
			break
		}
		err = os.Chmod(fileName, 0777)
		if err != nil {
			passed = false
			break
		}
	}

	// set and test env variable, then add to .profile
	os.Setenv("gopherSeq_bin", newBin)
	if len(os.Getenv("gopherSeq_bin")) == 0 {
		return false
	} else {
		exportCmd := "export gopherSeq_bin=\"" + newBin + "\"\n"
		f, err := os.OpenFile(homeDir+"/.profile", os.O_APPEND|os.O_WRONLY, 0600)
		if err != nil {
			panic(err)
		}
		defer f.Close()
		if _, err = f.WriteString(exportCmd); err != nil {
			fmt.Println("couldn't add gopherSeq_bin export statement to .profile file!\n\n")
			os.Exit(1)
		}
		return true
	}
}

/*
  functions to select which envtest to run
*/
func Test4align_progs() (bool, []string) {
	passed, messages = true, nil
	passed, messages = ProgramTest(align_programs)
	return passed, messages
}
func Test4qcheck_progs() (bool, []string) {
	passed, messages = true, nil
	passed, messages = ProgramTest(check_programs)
	return passed, messages
}

/*
  functions to test for installed software
*/
func ProgramTest(required_programs []string) (bool, []string) {
	for _, program := range required_programs {
		program_path, err := exec.Command("which", program).Output()
		if err != nil {
			messages = append(messages, " * can't find "+program+"!\n")
			passed = false
		} else {
			messages = append(messages, " * found "+program+" --> "+string(program_path))
		}

		// make sure they are using samtools version 1.4
		if program == "samtools" {
			program_path, err := exec.Command("samtools", "--version").Output()
			if err != nil {
				messages = append(messages, " * your version of samtools seems incorrect - please check for >= 1.4\n")
				passed = false
			} else if match, _ := regexp.MatchString("samtools 1.4", string(program_path)); match == false {
				messages = append(messages, " * please update your samtools to version 1.4\n")
				passed = false
			}
		}

		// make sure they are using bcftools version 1.4
		if program == "bcftools" {
			program_path, err := exec.Command("bcftools", "--version").Output()
			if err != nil {
				messages = append(messages, " * your version of bcftools seems incorrect - please check for >= 1.4\n")
				passed = false
			} else if match, _ := regexp.MatchString("bcftools 1.4", string(program_path)); match == false {
				messages = append(messages, " * please update your bcftools to version 1.4\n")
				passed = false
			}
		}
	}
	return passed, messages
}

///////////////
// MAIN
//////////////
func Main() {
	// print usage or parse arguments
	if len(os.Args) < 2 {
		printInfo()
	} else {
		arg.MustParse(&args)
	}

	// check for installed software
	if args.Run != false {
		fmt.Printf("testing for gopherSeq bin . . .\n")
		passed, messages = BinCheck()
		for _, message := range messages {
			fmt.Printf("%v", message)
		}
		if passed == false {
			fmt.Printf("\n!\nplease address the above gopherSeq_bin error!\n!\n\n")
			os.Exit(1)
		}
		fmt.Printf("testing for required software . . .\n")
		fmt.Printf("QC check programs:\n")
		passed, messages = Test4qcheck_progs()
		for _, message := range messages {
			fmt.Printf("%v", message)
		}
		if passed == false {
			fmt.Printf("\n!\nprogram test failed for required QC check programs!\n!\n\n")
		}
		fmt.Printf("Align programs:\n")
		passed, messages = Test4align_progs()
		for _, message := range messages {
			fmt.Printf("%v", message)
		}
		if passed == false {
			fmt.Printf("\n!\nprogram test failed for required align programs!\n!\n\n")
		}
	}
}
