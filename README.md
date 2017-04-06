[![Build Status](https://travis-ci.org/will-rowe/gopherSeq.svg?branch=master](https://travis-ci.com/will-rowe/gopherSeq)
<img src="https://github.com/will-rowe/will-rowe.github.io/raw/master/images/gopher.png" align="right" width="68" >

# gopherSeq

A set of simple pipelines and tools for microbial whole genome sequence data

***

`this is still under development - improvements and new tools are being added`

***

##  Installation

### Download a binary

Download, uncompress and add the binary to your PATH:

* [linux](https://github.com/will-rowe/gopherSeq/releases/download/0.0.1/gopherSeq.linux.tar)
* [OSX](https://github.com/will-rowe/gopherSeq/releases/download/0.0.1/gopherSeq.osx.tar)

### Alternatively

Install using Go:
```
go get -u github.com/will-rowe/gopherSeq/cmd/gopherSeq
go install github.com/will-rowe/gopherSeq/cmd/gopherSeq
```

Or, download the source and compile into executable:
```
git clone https://github.com/will-rowe/gopherSeq
cd ./gopherSeq/cmd/gopherSeq
env GOOS=linux GOARCH=amd64 go build   ## to compile for cluster
env GOOS=darwin GOARCH=amd64 go build  ## to compile for OSX
```

### Requirements

There are several pieces of software and environment variables needed before running this program. The program will check for environment variables and software requirements before starting any of the commands and it will let you know if any aren't installed properly.

This is a list of required software (not included in this repo):

| qcheck | align |
| ------------- | ------------- |
| java | bowtie2 |
| fastqc | samtools (1.4) |
| trimmomatic | bcftools(1.4) |
| kraken | |
| multiqc | |



**IMPORTANT** --> In addition to the above software, this program also requires a special bin to be set (called `gopherSeq_bin`). The program can set this up for you, just run the `envtest` command. Alternatively, download the bin from this repo and create an environment variable to point to it:

```
echo export gopherSeq_bin=\"/path/to/gopherSeq/bin\" >> ~/.profile
```

If you want the `qcheck` command to run Kraken and Trimmomatic adapter removal, you need to create 2 symbolic links in the gopherSeq_bin. The naming of these links is important:
```
ln -s /path/to/kraken/minikraken_20141208 $gopherSeq_bin/kraken_db
ln -s /path/to/adapters/TruSeq3-SE.fa $gopherSeq_bin/adapters.fa
```

### Caveats

* this program has been designed to work well with our typical *salmonella* WGS data. As a result, a lot of the SNP calling options (samtools, GATK, bcftools etc.) have been hard-coded (for now...) - if you want to adjust the parameters, you'll have to edit the code and recompile (same goes for the QC options)

* only basic file checking is carried out - mainly extension checks and those done by the called programs

* when grouping read files into samples, only the extensions `*_1.fastq` and `*_2.fastq` are considered (could also be .fq and/or gzipped). This means that singleton read files are not merged into samples


***

## Commands


### envtest

This just tests the environment to make sure that the required software can be found and that the `gopherSeq_bin` has been set up. If it needs to set up the bin, it will download the files, add an export statement to your .profile file and then exit.

Basic usage:
```
gopherSeq envtest --run
```

### qcheck

A *very* basic quality checking pipeline. This won't inspect any of the QC results, it just runs a series of QC programs and makes a pretty report with multiqc. You can go submit the trimmed reads straight to the `align` tool (using options --align and --reference ./xxx.fa). There is no log with this tool - all output is straight to STDOUT.

Basic usage:
```
gopherSeq qcheck /path/to/input/*.fastq.gz
```

### align

This is a simple pipeline for aligning and variant calling bacterial WGS data against a reference. It takes fastq reads and a reference, performs an alignment for each sample, runs GATK indel correction, calls SNPs and then creates a pseudogenome for each sample (for use in downstream phylogenetic analyses). This command can accept a mix of paired end data and single-end --> it stores paired-end data under a single sample name ONLY if the files end in `_1.fastq` (or variant e.g. `_1.fq.gz`)

Basic usage:
```
gopherSeq align --reference /path/to/reference.fasta /path/to/input/*.fastq.gz
```
