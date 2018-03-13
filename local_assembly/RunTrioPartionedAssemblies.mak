# configuration for paths, etc
include $(PBS)/local_assembly/Configure.mak

COVERAGE?=10

TMPNAME=$(shell echo $(REGION) | sed 's/:/_/')
STREGION=$(shell echo $(REGION) | sed 's/\./:/')

all: h0.sam h0/h0.fasta h1.sam h1/h1.fasta assembly.consensus.fasta.sam \
  assemblies samfiles results

help:
	echo "Usage: make -f RunDiploidAssembly.py REGION=chr:start-end"
	echo "REF=/path/to/reference.fasta "
	echo "BAMS=aligned_bams.fofn "
  echo "SAMASSEMBLER=/path/to/sam_assembler.mak "
  echo " CH=child-id "
  echo " MAKE_CONFIG=configfile.mak "
  echo " MOBAMS=path-to-mother-bams MO=mother-id MO_INHERIT=mo.inherit.bam.gz "
  echo " FABAMS=path-to-mother-bams fA=mother-id FA_INHERIT=mo.inherit.bam.gz "


region.vcf:
	tabix -h $(VCF) $(STREGION) > $@

h0.sam: region.vcf
	head -1 $(CHBAMS) | xargs -i samtools view -H {} > reads.sam
	cat $(CHBAMS) | xargs -i samtools view -q 30 {} $(STREGION)  |  sed -e "s/qi/iq/" -e "s/qd/dq/" -e "s/qs/sq/" -e "s/qm/mq/" -e "s/td/dt/" -e "s/ts/st/" >> reads.sam
	$(PBS)/local_assembly/samToBed reads.sam | $(PBS)/local_assembly/DetectChimeras.py > filter.list
	grep -v -f filter.list reads.sam | $(PBS)/local_assembly/RemoveShortSubreads.py 500 | $(PBG)/partitionByPhasedSNVs --vcf region.vcf --sam /dev/stdin --rgn $(REGION) --pad 100000 --h1 h0.sam --h2 h1.sam --ref $(REF) --minGenotyped 1 --summary summary.txt --sample $(CH) --unassigned unassigned.sam --minScoreDifference 1 --nw-window 7
	-cp h0.sam h0.sam.bak
	-cp h1.sam h1.sam.bak
	$(PBS)/local_assembly/ConcatenateUnassignedReads.sh $(PBS)/local_assembly/ConcatenateUnassignedReads.py $(VCF) $(CH) $(REGION) 20 unassigned.sam  h0.sam h1.sam

  # partition parents now
	echo "Selecting from mother"
	$(PBS)/local_assembly/SelectInherited.sh \
   -b $(MOBAMS) \
   -v $(VCF) \
   -i $(MO_INHERIT) \
   -r $(STREGION) \
   -m mo \
   -R $(REF) \
   -s $(MO)

	-grep -v "^@" mo.inherited.sam >> h`cat mo.chrom.txt|tr -d "\n"`.sam
	echo "catted  mo.inherited.sam to h`cat mo.chrom.txt`.sam "
	echo "Selecting from father."
  # partition parents now
	$(PBS)/local_assembly/SelectInherited.sh \
   -b $(FABAMS) \
   -v $(VCF) \
   -i $(FA_INHERIT) \
   -r $(STREGION) \
   -m fa \
   -R $(REF) \
   -s $(FA)
	-grep -v "^@" fa.inherited.sam >> h`cat fa.chrom.txt | awk '{print $2;}'|tr -d "\n"`.sam



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
	cd h0 && source $(PBS)/local_assembly/setup_assembly.sh && make -f $(SAMASSEMBLER) SAM=../h0.sam HAP=H0 REGION=$(STREGION) REF=$(REF) COV=$(COVERAGE) && mv assembly.consensus.fasta h0.fasta  && mv assembly.consensus.fasta.sam ../h0.var.sam 

h1/h1.fasta: h1.sam
	mkdir -p h1
	cd h1 && source $(PBS)/local_assembly/setup_assembly.sh && make -f $(SAMASSEMBLER) SAM=../h1.sam HAP=H1 REGION=$(STREGION) REF=$(REF) COV=$(COVERAGE) && mv assembly.consensus.fasta h1.fasta && mv assembly.consensus.fasta.sam ../h1.var.sam 

assembly.consensus.fasta.sam: h0.var.sam h1.var.sam
	cat h0.var.sam | awk '{ print $$0"\tHA:i:1";}' > $@
	grep -v "^@" h1.var.sam | awk '{ print $$0"\tHA:i:2";}' >> $@
	cat h0/h0.fasta | tr "|" "_" | awk '{ if (substr($$1,0,1) == ">") {print $$1"/1";} else{ print $$1;} }' > assembly.consensus.fasta
	cat h1/h1.fasta | tr "|" "_" | awk '{ if (substr($$1,0,1) == ">") {print $$1"/2";} else{ print $$1;} }' >> assembly.consensus.fasta

clean:
	echo "not cleaning"
	#rm -rf 0-rawreads/    2-asm-falcon/	 alignment.sam		   assembly.fasta      input.fofn  reads.bas.h5  scripts/ 1-preads_ovl/  alignment.cmp.h5  assembly.fasta.fai  reads.bam   reads.fasta	 sge_log/ assembly.bam aligend_reads.vcf aligned_reads.bam aligned_reads.bam.bai aligned_reads.bam.nft  assembly.consensus.fasta.fai  h1 h2 reads.sam 

