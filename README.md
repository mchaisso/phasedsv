Phased-SV
=========

Local assembly based SV calling with single-molecule sequencing reads
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

1. Quiver.

Install quiver binary scripts into quiver/bin, and supporting
libraries into quiver/lib.  Because this can be a cumbersome task
depending on your software environment, precompiled binaries are
distributed in quiver_bin. If you wish to avoid trying to install
quiver/arrow from source, move the quiver_bin directory into quiver:

```mv quiver_bin quiver```





