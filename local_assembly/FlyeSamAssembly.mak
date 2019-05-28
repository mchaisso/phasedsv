TMPNAME=$(shell echo $(REGION) | sed 's/:/_/')
MIN_QUALITY?=20
MAKE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

SHELL=/bin/bash

# other setup
CWD=$(shell pwd)

all: reads.fasta assembly.consensus.fasta assembly.consensus.fasta.sam assembly.fasta assembly.bam

test:
	echo "$(inOptions)"

help:
	@echo "Usage make -f FlyeSamAssembly.mak SAM=samfile [REF=reference.fasta] "


reads.fasta: $(SAM)
	grep -v "^@" $(SAM) | $(MAKE_DIR)/StreamSamToFasta.py | $(MAKE_DIR)/FormatFasta.py --fakename  > $@


assembly.fasta: reads.fasta
	. $(MAKE_DIR)/../dep/build/bin/activate python2 && flye --pacbio-raw reads.fasta --genome-size 70000 -o flye -t 4 -i 0
	cp flye/assembly.fasta $@

assembly.bam: assembly.fasta $(SAM)
	export READ_SOURCE=$(READ_SOURCE) && $(MAKE_DIR)/MakeBamVersionWrapper.sh $(SAM) $(READ_SOURCE) assembly.fasta assembly.bam
	samtools index assembly.bam


assembly.bam.pbi: assembly.bam
	pbindex assembly.bam

assembly.consensus.fasta: assembly.bam  assembly.fasta
	samtools faidx assembly.fasta
	echo $(PYTHONPATH)
	$(MAKE_DIR)/RunConsensusTool.sh assembly.fasta assembly.bam $@ $(READ_SOURCE)


#assembly.consensus.fasta: assembly.bam assembly.bam.pbi assembly.fasta
#	samtools faidx assembly.fasta
#	echo $(PYTHONPATH)
#	$(MAKE_DIR)/RunConsensusTool.sh assembly.fasta assembly.bam $@ $(READ_SOURCE)
#	awk '{ if (substr($$1,0,1) == ">") {print $$1"/$(HAP)";} else { print;} }' $@ > $@.tmp
#	mv -f $@.tmp $@

assembly.consensus.fasta.sam: assembly.consensus.fasta
	echo $(REGION) | tr ":-" "\t\t" > region.bed
	bedtools slop -i region.bed -g $(REF).fai -b 10000 | awk '{ print $$1":"$$2"-"$$3;}' > region.wide
	region=`cat region.wide`
	samtools faidx $(REF) `cat region.wide` > target.ref.fasta
	$(MAKE_DIR)/../dep/bin/blasrmc assembly.consensus.fasta target.ref.fasta -sam -bestn 1 -maxMatch 30 -sdpTupleSize 9 -indelRate 3 -affineAlign -affineOpen 8 -affineExtend 0  -clipping soft | $(MAKE_DIR)/shiftSamPos > $@

clean:
	echo "Not cleaning up after Canu"
	#rm -rf 0-rawreads/ 1-preads_ovl   2-asm-falcon/	 alignment.sam		   assembly.fasta      input.fofn  reads.bas.h5  scripts/ 1-preads_ovl/  alignment.cmp.h5  assembly.fasta.fai  reads.bam   reads.fasta	 sge_log/


