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

Installation
-----------

git clone --recursive https://github.com/mchaisso/phasedsv.git


1. Required software. This should be default software accessible by your path.

    samtools >= 1.2 ( htslib >= 1.2.1)
		tabix >= 0.2.5
		bedtools >= 2.27.1
    vt >= 0.5772
		ucsc software (bedToBigBed)
    R
		virtualenv
		canu 

2. Install quiver.

Install quiver binary scripts into quiver/bin, and supporting
libraries into quiver/lib.  Because this can be a cumbersome task
depending on your software environment, precompiled binaries are
distributed in quiver_bin. If you wish to avoid trying to install
quiver/arrow from source, move the quiver_bin directory into quiver with the commands below. This is built for 64-bit CentOS.

```mv quiver_bin quiver
export PATH=$PATH:$PWD/quiver/bin
export PYTHONPATH=$PYTHONPATH:$PWD/quiver/lib/python2.7/site-packages/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/quiver/lib
```

3. Build binary requirements.

```make```

4. Add python requirements

```cd hgsvg && ./setup_virtualenv.sh```



