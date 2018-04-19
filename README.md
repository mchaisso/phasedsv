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
    2.1.  Setup python environment. It is easiest to ensure compatibility with python modules if virtual environments are used.  A utility script is provided to create and populate the virtualenv with the required python modules for running phasedsv.  You can create the module using `source phasedsv/setup_virtualenv.sh`

    2.2. `phasedsv/setup_phasedsv.sh` : configuration of the python environment,
  and the LD_LIBRARY_PATH or PATH to access samtools and bedtools.  A
  sample script is included.

    2.3 `phasedsv/local_assembly/Configure.mak`:  This sets up variables
  used in the make files that run local assemblies. They need to point
  to the reference (indexed by samtools faidx and blasr sawriter), and
  the canu installation. The value of "READ_SOURCE" should be set to
  "HGSVG_BAM" if running PacBio RSII alignments, and anything else (or
  not set) otherwise.

    2.4 BAM fofn: this should be a file of complete paths to the bam or
	bams if the alignments are split into multiple bams.

    2.4 Trio assembly.
	   If you are running a trio assembly, you need to generate bed
	   files that contain the regions which the parental reads may be
	   unambiguously assigned.  The example below uses the phased vcf
	   for the Puerto Rican family.

`phasedsv/hgsvg/phasing/DetermineInheritance.py  --vcf bams/PR05.wgs.whatshap.strandseq-10X.20160704.phased-genotypes.vcf.gz --child HG00733 --fa HG00731 --mo HG00732 --faBed fa.bed --moBed mo.bed`


3. Run local assemblies.

  3.1 Create a parameter file describing the source data for the file. This requires the following values
. The full path to the reference reads are aligned to.
. The path to the bam file of file names. This is one line per bam, with the full path to the file.
. The VCF file that has the phased SNVs. 
. The name of the sample that is being assembled (the sample ID in the VCF file). 

```
export REF=/home/cmb-panasas2/mchaisso/genomes/hg38/hg38.fa
export BAMS=/home/cmb-panasas2/mchaisso/giab/bams.fofn
export VCF=/home/cmb-16/mjc/sample-datasets/giab/HG002/pacbio/10x/NA24385_GRCh38.het.vcf.gz
export SAMPLE=HG002
```

 3.2 Define the regions that will be assembled. These can be copied from `phasedsv/hgsvg/regions/Windows.60kb-span.20kbp-stride.txt.  In general they are in the format chrom.start-end. 

 3.3 Run local assemblies. This is a computationally intensive task, and is best ran on a cluster. For every line in the regions file, a local assembly must be generated. For single-sample assemblies, the command is RunTiledAssembly.sh

```
Usage: $base/RunTiledAssembly.sh region assembly.params [options]

Options: -j STR   Name of job running on grid engine.
         -a FILE  Full path to local assembler file, for example 
                  $BASE/local_assembly/CanuSamAssembly.mak
         -d STR   Working directory.
```

For trio assemblies, you can run RunTrioTiledAssemblyOnRegions.sh
```
Usage $base/RunTrioTiledAssemblyOnRegions.sh region assembly.params [options]
Options -j STR    Name of job running on grid engine.
```

Because this is many small jobs, SGE and PBS work best with array jobs. There are sample jobs for PBS for single-genome at phasedsv 




