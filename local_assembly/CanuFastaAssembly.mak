# configuration for paths, etc
PBS=/home/cmb-16/mjc/mchaisso/projects/LocalAssembly/scripts

include $(PBS)/local_assembly/Configure.mak



#all: assembly.fasta reads.bam reads.fasta

MIN_QUALITY?=20



# other setup
CWD=$(shell pwd)

all: assembly.fasta
	@echo $(BTREGION)

test:
	echo "$(inOptions)"

help:
	@echo "Usage make -f CanuFastaAssembly.mak READS=file.fasta [REF=reference.fasta] [BAS=one_base_file]"

PBS=/home/cmb-16/mjc/mchaisso/projects/LocalAssembly/scripts

CANU_DIR=/home/cmb-16/mjc/shared/software_packages/canu/Linux-amd64/bin/

assembly.fasta: $(READS)
#	$(PBS)/FastaToFakeFastq.py $(READS) reads.fastq
	source /usr/usc/java/1.8.0_45/setup.sh && $(CANU_DIR)/canu -nanopore-raw $(READS) genomeSize=60000 -d assembly -p asm useGrid=false  gnuplotTested=true   cnsThreads=2 ovlThreads=2 mhapThreads=4 
	cp assembly/asm.contigs.fasta $@

	-rm -rf templocal


clean:
	rm -rf reads.fasta 0-rawreads/ 1-preads_ovl/ 2-asm-falcon/ scripts/ sge_log/
