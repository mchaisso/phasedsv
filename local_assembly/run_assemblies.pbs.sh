#!/usr/bin/env bash
#PBS -t 1-148923%150
#PBS -l walltime=30:00:00,nodes=1:ppn=4,mem=4gb
base=dirname $0

$base/RunTiledAssembly.sh `awk  "NR==$PBS_ARRAYID" < Windows.60kb-span.20kbp-stride.txt` assembly.config





