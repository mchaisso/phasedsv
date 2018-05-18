Installing Phased-SV.

1. Requirements:

  gcc >= 7.0
  CMAKE >= 1.10
	ninja >= 1.8

These should all be accessible in your path.

2. Make binaries. The full tree and dependencies will be built by typing

`make` in the phasedsv directory.  This takes some time, and necessary
modules will be downloaded. 

It may be necessary to configure your build. In particular, if the g++ >= 7.0 is not the default directory, set the environment variables:

CC=/path/to/gcc
CXX=/path/to/g++


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
  
   The default configfile that is build for phasedsv,
   `setup_phasedsv.sh`, is set up to use the distributed `quiver_bin`
   after it is symlinked to `quiver`. 

