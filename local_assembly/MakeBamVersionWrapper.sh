#!/usr/bin/env bash

BASE=$(dirname $0)
SAM=$1
if [ $READ_SOURCE = "HGSVG_BAM" ]; then
	echo "USING HGSVG BAM"
		$BASE/blasr/pbihdfutils/bin/samtobas $SAM $SAM.bas.h5
		$BASE/blasr/alignment/bin/blasr $SAM.bas.h5 assembly.fasta  \
				-clipping subread -sam -nproc 8 -out  /dev/stdout -preserveReadTitle | \
				samtools view -bS - | \
				samtools sort -T tmp -o assembly.bam
else
	$BASE/blasr/alignment/bin/blasr $SAM assembly.fasta  \
        -fileType sam -clipping subread -passthrough -sam -nproc 8 -out  /dev/stdout -preserveReadTitle | \
        $BASE/FixReadGroup.py | \
        $BASE/../pbsamstream/pbsamstream  - | \
        $BASE/../samtools/samtools view -bS - | $BASE/../samtools/samtools sort -T tmp -o assembly.bam
fi

		

