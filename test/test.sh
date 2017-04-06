#!/bin/bash
go build -o ./gopherSeq ../cmd/gopherSeq/main.go

export gopherSeq_bin=../bin

./gopherSeq version
./gopherSeq qcheck ./data/reads/ERR1107833_downsampled_singletons.fastq.gz
./gopherSeq align --reference ./data/RefSeq/NC_004741.fasta ./data/reads/ERR1107833_downsampled_pass*.fastq.gz
