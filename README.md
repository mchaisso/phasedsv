Phased-SV
=========

Local assembly based SV detection using single-molecule sequencing reads
and a phased SNV VCF file.

Summary
-------

This software pipeline performs SV calling with four main steps:
1. Local assembly of haplotype-partitioned reads.
2. Merging of local assemblies into a reference-guided assembly.
3. Mapping merged assemblies to the reference.
4. Filtering SVs by read-support of breakpoints.


For installation, please refer to Install.md


Running PhasedSV
----------------

1. Data generation.

Phased-SV assumes you have BAM files of reads aligned to the
reference, and a vcf file of phased SNVs. A minimum of 40X works
best. To run a test on chromosome 22, you can download data listed in
SampleData.txt. If you are generating bams from scratch, you can
reference the README.md in pbsamstream for the commands to generate
correctly formatted bam files.

2. Configuration.

  2.1. `phasedsv/config.sh` : configuration of the python environment,
  and the LD_LIBRARY_PATH or PATH to access samtools and bedtools.  A
  sample script is included.

  2.2 `phasedsv/local_assembly/Configure.mak`:  This sets up variables
  used in the make files that run local assemblies. They need to point
  to the reference (indexed by samtools faidx and blasr sawriter), and
  the canu installation.

	2.3 BAM fofn: this should be a file of complete paths to the bam or
	bams if the alignments are split into multiple bams.

  2.4 Trio assembly.
	   If you are running a trio assembly, you need to generate bed
	   files that contain the regions which the parental reads may be
	   unambiguously assigned.  The example below uses the phased vcf
	   for the Puerto Rican family.

`~/projects/phasedsv/hgsvg/phasing/DetermineInheritance.py  --vcf bams/PR05.wgs.whatshap.strandseq-10X.20160704.phased-genotypes.vcf.gz --child HG00733 --fa HG00731 --mo HG00732 --faBed fa.bed --moBed mo.bed`


3. Run local assemblies.


