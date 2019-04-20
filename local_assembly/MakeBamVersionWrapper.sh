#!/usr/bin/env bash

BASE=$(dirname $0)
SAM=$1
READ_SOURCE=$2
pbmm2 index assembly.fasta assembly.fasta.mmi
a=$(basename $SAM)
b=${a%.*}

if [ "$READ_SOURCE"="HGSVG_BAM" ]; then
    echo "cat $SAM | $BASE/ReformatHGSVGSam.py | samtools view -bhS - -o $b.bam"
    cat $SAM | $BASE/ReformatHGSVGSam.py | samtools view -bhS - -o $b.bam
else
    samtools view -bhS $SAM -o $b.bam
fi
pbmm2 align --sort assembly.fasta.mmi $b.bam assembly.bam 
pbindex assembly.bam


		

