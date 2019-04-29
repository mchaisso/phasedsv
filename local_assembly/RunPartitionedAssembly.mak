# configuration for paths, etc
MAKE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))



COVERAGE?=10
CHROM=$(shell echo $(REGION) | tr ":-" "\t\t" | cut -f 1)
START=$(shell echo $(REGION) | tr ":-" "\t\t" | cut -f 2)
END=$(shell echo $(REGION) | tr ":-" "\t\t" | cut -f 3)
SNVSTART=$($START+$WINDOW)
SNVEND=$($END-$WINDOW)
SAMASSEMBLER=$(MAKE_DIR)/CanuSamAssembly.mak

TMPNAME=$(shell echo $(REGION) | sed 's/:/_/')
STREGION=$(shell echo $(REGION) | sed 's/\./:/')

all:  region.vcf h1.sam h2/assembly.consensus.fasta.sam h2/assembly.consensus.fasta.sam  h1.var.sam assembly.consensus.fasta.sam

help:
	echo "Usage: make -f RunPartitionedAssembly.mak " ;\
    echo "               REGION=chr:start-end "; \
    echo "               REF=/path/to/reference.fasta " ;\
    echo "               BAMS=aligned_bams.fofn" ;\
    echo "               VCF=/path/to/phased.vcf" ; \
    echo "               SAMASSEMBLER=/path/to/sam_assembler.mak " \
    echo "                          " ; \
	  echo "   optional parameters:" ;\
	  echo "   COVERAGE=cov    Minimum coverage to allow het call.";\
	  echo "   SUBASSEMBLER=/path/to/asm.mak   Assembler to call on sam files." 

region.vcf:
	tabix -h $(VCF) $(STREGION) > $@

h1.sam: region.vcf reads.sam 
	$(MAKE_DIR)/pbgreedyphase/partitionByPhasedSNVs --vcf region.vcf --sam reads.sam --h1 h1.sam --h2 h2.sam --rgn $(REGION) --ref $(REF) --nw-window 5 --minGenotyped 1 $(AUTO) --sample $(SAMPLE) >& summary.txt
	echo "h1" >> summary.txt
	grep  "^@" h1.sam | wc -l >> summary.txt
	echo "h2" >> summary.txt
	grep  "^@" h2.sam | wc -l >> summary.txt


reads.sam:
	samtools view -H `head -1 $(BAMS)` > $@.full
	cat $(BAMS) | xargs -i samtools view  -q 30 {} $(CHROM):$(START)-$(END) >> $@.full
	samToBed $@.full | $(MAKE_DIR)/DetectChimeras.py > filter.list
	grep -v -f filter.list $@.full > $@


h1/assembly.consensus.fasta.sam: h1.sam
	mkdir -p h1
	cd h1 && make -f $(SAMASSEMBLER) SAM=../h1.sam HAP=H1 REGION=$(STREGION) REF=$(REF) || true

h2/assembly.consensus.fasta.sam: h2.sam h1/assembly.consensus.fasta.sam
	mkdir -p h2
	if [ $(IS_AUTO) == "AUTO" ]; then \
    cp h1/assembly.consensus.fasta h2/; \
    cp h1/assembly.consensus.fasta.sam h2/; \
  else \
	 cd h2 && make -f $(SAMASSEMBLER) SAM=../h2.sam HAP=H2 REGION=$(STREGION) REF=$(REF) || true ; \
  fi

h1.var.sam: h1/assembly.consensus.fasta.sam
	touch $@
	-cat h1/assembly.consensus.fasta.sam | awk '{ $$1=$$1"/1"; print; }' | tr " " "\t" > h1.var.sam 
	-cp h1/assembly.consensus.fasta h1.fasta

h2.var.sam: h2/assembly.consensus.fasta.sam
	touch $@
	-cat h2/assembly.consensus.fasta.sam  | awk '{ $$1=$$1"/2"; print; }' | tr " " "\t" > h2.var.sam 
	-cp h2/assembly.consensus.fasta h2.fasta


assembly.consensus.fasta.sam: h1.var.sam h2.var.sam
	cp h1.var.sam $@
	grep -v "^@" h2.var.sam >> $@
	-cat h1.fasta | tr "|" "_" | awk '{ if (substr($$1,0,1) == ">") {print $$1"/1";} else{ print $$1;} }' > assembly.consensus.fasta
	-cat h2.fasta | tr "|" "_" | awk '{ if (substr($$1,0,1) == ">") {print $$1"/2";} else{ print $$1;} }' >> assembly.consensus.fasta

clean:
	rm -rf 0-rawreads/    2-asm-falcon/	 alignment.sam		   assembly.fasta      input.fofn  reads.bas.h5  scripts/ 1-preads_ovl/  alignment.cmp.h5  assembly.fasta.fai  reads.bam   reads.fasta	 sge_log/ assembly.bam aligend_reads.vcf aligned_reads.bam aligned_reads.bam.bai aligned_reads.bam.nft  assembly.consensus.fasta.fai  h1 h2 reads.sam 

