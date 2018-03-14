include Configure.mak
TMPNAME=$(shell echo $(REGION) | sed 's/:/_/')
MIN_QUALITY?=20
MAKE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# other setup
CWD=$(shell pwd)

all: reads.fasta assembly.fasta assembly.bam assembly.consensus.fasta assembly.consensus.fasta.sam

test:
	echo "$(inOptions)"

help:
	@echo "Usage make -f CanuSamAssembly.mak SAM=samfile [REF=reference.fasta] "


reads.fasta: $(SAM)
	grep -v "^@" $(SAM) | $(PBS)/local_assembly/StreamSamToFasta.py | $(PBS)/local_assembly/FormatFasta.py --fakename  > $@

# corMhapSensitivity=high corMinCoverage=2 errorRate=0.035 trimReadsCoverage=4 cnsMemory=2G cnsThreads=1  corMhapSensitivity=high corMinCoverage=3 errorRate=0.025 errorRate=0.10 utgOverlapper=ovl utgOvlErrorRate=0.15
assembly.fasta: reads.fasta
	$(CANU_DIR)/canu -pacbio-raw reads.fasta genomeSize=60000 -d assembly -p asm useGrid=false  gnuplotTested=true  corMhapSensitivity=high corMinCoverage=1 cnsThreads=4 ovlThreads=4 mhapThreads=4 contigFilter="2 1000 1.0 1.0 2"
	if [ -s assembly/asm.contigs.fasta ]; then \
    cp assembly/asm.contigs.fasta $@; \
  fi
#	if [ ! -s assembly/asm.contigs.fasta ]; then \
#		ls reads.fasta  > input.fofn; \
#		source $(HOME)/scripts/setup_falcon.sh && fc_run.py $(HOME)/projects/PacBioSequencing/scripts/local_assembly/falcon/fc_run.local.cfg; cp 2-asm-falcon/p_ctg.fa $@;   \
#	else \
#	   cp assembly/asm.contigs.fasta $@; \
#        fi ;
#	echo "Number of reads " > report.txt
#	grep -c ">" reads.fasta >> report.txt
#	echo "Assembly number of contigs" >> report.txt
#	module load numpy/latest; ~mchaisso/software/mcsrc/UTILS/pcl reads.fasta | ~/scripts/stats.py >> report.txt
#	-rm -rf templocal

assembly.bam: assembly.fasta 
	$(MAKE_DIR)/../hgsvg/alignment/bin/blasr $(SAM) assembly.fasta  \
        -clipping subread -passthrough -sam -nproc 8 -out  /dev/stdout -preserveReadTitle | \
        $(PBS)/pbsamstream/pbsamstream  - | \
         $(PBS)/samtools/samtools view -bS - | $(PBS)/samtools/samtools sort -T tmp -o assembly.bam
	samtools index assembly.bam

assembly.bam.pbi: assembly.bam
	source $(DEPLOY)/setup-env.sh && $(DEPLOY)/bin/pbindex assembly.bam

assembly.consensus.fasta: assembly.bam assembly.bam.pbi assembly.fasta
	samtools faidx assembly.fasta
	source $(DEPLOY)/setup-env.sh && export LD_LIBRARY_PATH=/usr/usc/gnu/gcc/5.3.0/lib64/:$$LD_LIBRARY_PATH && $(DEPLOY)/bin/quiver  -j4 --minCoverage 7 --noEvidenceConsensusCall nocall --referenceFilename assembly.fasta assembly.bam -o $@ 
	awk '{ if (substr($$1,0,1) == ">") {print $$1"/$(HAP)";} else { print;} }' $@ > $@.tmp
	mv -f $@.tmp $@

assembly.consensus.fasta.sam: assembly.consensus.fasta
	echo $(REGION) | tr ":-" "\t\t" > region.bed
	bedtools slop -i region.bed -g $(REF).fai -b 10000 | awk '{ print $$1":"$$2"-"$$3;}' > region.wide
	region=`cat region.wide`
	samtools faidx $(REF) `cat region.wide` > target.ref.fasta
	$(MAKE_DIR/)/alignment/bin/blasr assembly.consensus.fasta target.ref.fasta -sam -bestn 1 -maxMatch 30 -sdpTupleSize 9 -indelRate 3 -affineAlign -affineOpen 8 -affineExtend 0  -clipping soft | $(PBS)/local_assembly/shiftSamPos > $@

clean:
	echo "Not cleaning up after Canu"
	#rm -rf 0-rawreads/ 1-preads_ovl   2-asm-falcon/	 alignment.sam		   assembly.fasta      input.fofn  reads.bas.h5  scripts/ 1-preads_ovl/  alignment.cmp.h5  assembly.fasta.fai  reads.bam   reads.fasta	 sge_log/


