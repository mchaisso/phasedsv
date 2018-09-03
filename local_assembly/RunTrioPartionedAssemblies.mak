# configuration for paths, etc

MAKE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(MAKE_DIR)/Configure.mak


COVERAGE?=10
SHELL=/bin/bash
TMPNAME=$(shell echo $(REGION) | sed 's/:/_/')
STREGION=$(shell echo $(REGION) | sed 's/\./:/')
SAMASSEMBLER=$(MAKE_DIR)/CanuSamAssembly.mak

all: h0.sam h0/h0.fasta h1.sam h1/h1.fasta assembly.consensus.fasta.sam \
  assemblies samfiles results

help:
	echo "Usage: make -f RunDiploidAssembly.py REGION=chr:start-end"
	echo "REF=/path/to/reference.fasta "
	echo "BAMS=aligned_bams.fofn "
	echo "SAMASSEMBLER=/path/to/sam_assembler.mak "
	echo " SAMPLE=child-id "
	echo " MAKE_CONFIG=configfile.mak "
	echo " MOBAMS=path-to-mother-bams MO=mother-id MO_INHERIT=mo.inherit.bam.gz "
	echo " FABAMS=path-to-mother-bams fA=mother-id FA_INHERIT=mo.inherit.bam.gz "


region.vcf:
	tabix -h $(VCF) $(STREGION) > $@

h0.sam: region.vcf
	head -1 $(BAMS) | xargs -i samtools view -H {} > reads.sam
	cat $(BAMS) | xargs -i samtools view -q 30 {} $(STREGION)  |  sed -e "s/qi/iq/" -e "s/qd/dq/" -e "s/qs/sq/" -e "s/qm/mq/" -e "s/td/dt/" -e "s/ts/st/" >> reads.sam
	if [ "$(IS_AUTO)" == "AUTO" ]; then \
    cp reads.sam h0.sam ; \
    cp reads.sam h1.sam ; \
  else \
    $(MAKE_DIR)/../mcutils/src/samToBed reads.sam | $(MAKE_DIR)/DetectChimeras.py > filter.list;  \
    grep -v -f filter.list reads.sam | $(MAKE_DIR)/RemoveShortSubreads.py 500 | $(MAKE_DIR)/pbgreedyphase/partitionByPhasedSNVs --vcf region.vcf --sam /dev/stdin --rgn $(REGION) --pad 100000 --h1 h0.sam --h2 h1.sam --ref $(REF) --minGenotyped 1 --summary summary.txt --sample $(SAMPLE) --unassigned unassigned.sam --minScoreDifference 1 --nw-window 7; \
  $(MAKE_DIR)/ConcatenateUnassignedReads.sh $(MAKE_DIR)/ConcatenateUnassignedReads.py $(VCF) $(SAMPLE) $(REGION) 20 unassigned.sam  h0.sam h1.sam; \
  fi 
  # partition parents now
  $(MAKE_DIR)/SelectInherited.sh \
   -b $(MOBAMS) \
   -v $(VCF) \
   -i $(MO_INHERIT) \
   -r $(STREGION) \
   -m mo \
   -R $(REF) \
   -s $(MO) ; \
  grep -v "^@" mo.inherited.sam >> h`cat mo.chrom.txt|tr -d "\n"`.sam ; \
  $(MAKE_DIR)/SelectInherited.sh \
   -b $(FABAMS) \
   -v $(VCF) \
   -i $(FA_INHERIT) \
   -r $(STREGION) \
   -m fa \
   -R $(REF) \
   -s $(FA); \
  grep -v "^@" fa.inherited.sam >> h`cat fa.chrom.txt | awk '{print $2;}'|tr -d "\n"`.sam


assemblies:
	mkdir -p assemblies

samfiles:
	mkdir -p samfiles

results:
	mkdir -p results

# h1 and h0 are generated at the same time, so just update it
h1.sam: h0.sam
	touch h1.sam

h0/h0.fasta: h0.sam
	mkdir -p h0
	cd h0 && source $(MAKE_DIR)/setup_assembly.sh && make -f $(SAMASSEMBLER) SAM=../h0.sam HAP=H0 REGION=$(STREGION) REF=$(REF) COV=$(COVERAGE) && mv assembly.consensus.fasta h0.fasta  && mv assembly.consensus.fasta.sam ../h0.var.sam 

h1/h1.fasta: h1.sam h0/h0.fasta
	mkdir -p h1;
	if [ $(IS_AUTO) == "AUTO" ]; then  \
		cp h0/h0.fasta h1/h1.fasta; \
  else \
    cd h1 && source $(MAKE_DIR)/setup_assembly.sh && make -f $(SAMASSEMBLER) SAM=../h1.sam HAP=H1 REGION=$(STREGION) REF=$(REF) COV=$(COVERAGE) && mv assembly.consensus.fasta h1.fasta && mv assembly.consensus.fasta.sam ../h1.var.sam; \
  fi


assembly.consensus.fasta.sam: h0.var.sam h1.var.sam
	cat h0.var.sam | awk '{ print $$0"\tHA:i:1";}' > $@
	grep -v "^@" h1.var.sam | awk '{ print $$0"\tHA:i:2";}' >> $@
	cat h0/h0.fasta | tr "|" "_" | awk '{ if (substr($$1,0,1) == ">") {print $$1"/1";} else{ print $$1;} }' > assembly.consensus.fasta
	cat h1/h1.fasta | tr "|" "_" | awk '{ if (substr($$1,0,1) == ">") {print $$1"/2";} else{ print $$1;} }' >> assembly.consensus.fasta

clean:
	rm -rf 0-rawreads/    2-asm-falcon/	 alignment.sam		   assembly.fasta      input.fofn  reads.bas.h5  scripts/ 1-preads_ovl/  alignment.cmp.h5  assembly.fasta.fai  reads.bam   reads.fasta	 sge_log/ assembly.bam aligend_reads.vcf aligned_reads.bam aligned_reads.bam.bai aligned_reads.bam.nft  assembly.consensus.fasta.fai  h1 h2 reads.sam 

