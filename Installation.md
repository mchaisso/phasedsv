Installing Phased-SV.

1. Requirements:

  gcc >= 4.9.0
  CMAKE >= 1.6
	ninja >= 1.8

These should all be accessible in your path.

2. Make binaries. The full tree and dependencies will be built by typing

`make` in the phasedsv directory.  This takes some time, and necessary
modules will be downloaded. 

3. Prepare quiver.
   An executable for quiver must exist in quiver/bin/quver. You can
   build this by cloning pacbio's pitchfork, and building
   GenomicConsensus, and setting up the directory structure, but this
   is often a difficult task.  A binary distribution of quiver is
   included, though it may not work on all systems. To try out the
   binary installation, link quiver_bin into quiver `ln -s quiver_bin
   quiver`, and then source setup_phasedsv.sh. Next, try typing
   `quiver` to see if paths are correctly configured, and there are no
   conflicts with versions of python, etc.

