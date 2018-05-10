all: shiftSamPos pbgreedyphase/partitionByPhasedSNVs blasr/alignment/bin/blasr

MAKE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

shiftSamPos: ShiftSamFromRegion.cpp
	g++ $< -o $@

blasr/alignment/bin/blasr:
	cd blasr && make HDF5INCLUDEDIR=$(MAKE_DIR)/../hdf5/build/include  HDF5LIBDIR=$(MAKE_DIR)/../hdf5/build/lib -j 4

mcutils/src/samToBed:
	cd mcutils/src && make

pbgreedyphase/partitionByPhasedSNVs:
	cd pbgreedyphase && make
