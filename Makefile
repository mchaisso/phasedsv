all: local_assembly/shiftSamPos \
	mcutils/src/samToBed \
  local_assembly/pbgreedyphase/partitionByPhasedSNVs \
  local_assembly/blasr/alignment/bin/blasr \
  hdf5/build/lib/libhdf5_cpp.a \
  hgsvg/blasr/alignment/bin/blasr \
  hgsvg/blasr/pbihdfutils/bin/samtobas \
  pbsamstream/pbsamstream \
  samtools/samtools \
  environments/python2.7/bin/activate \
  setup_phasedsv.sh

environments/python2.7/bin/activate:
	./setup_virtualenv.sh

setup_phasedsv.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "source "$(PWD)"/environments/python2.7/bin/activate" >> $@
	echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:"$(PWD)"/local_assembly/pbgreedyphase/boost_1_66_0/stage/lib/">> $@
	echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:"$(PWD)"/quiver/lib/">> $@
	echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:"$(PWD)"/hdf5/build/lib/">> $@
	echo "export PYTHONPATH=\$$PYTHONPATH:"$(PWD)"/quiver/lib/python2.7/site-packages/">> $@
	echo "export PATH=\$$PATH:"$(PWD)"/quiver/bin/" >> $@
	echo "#" >> $@
	echo "# Add custom configuration here." >> $@
	echo "#" >> $@

local_assembly/shiftSamPos:
	cd local_assembly && make

mcutils/src/samToBed:
	cd mcutils/src && make -j 8

local_assembly/pbgreedyphase/partitionByPhasedSNVs:
	cd local_assembly && make

local_assembly/blasr/alignment/bin/blasr: hdf5/build/lib/libhdf5_cpp.a
	cd local_assembly && make


hdf5/build/lib/libhdf5_cpp.a:
	mkdir -p $(PWD)/hdf5/build
	cd hdf5/ && \
  mkdir cmake_build && \
  cd cmake_build && \
  CXXFLAGS=-std=c++11 && cmake .. -DCMAKE_C_COMPILER=`which gcc` -DCMAKE_CPP_COMPILER=`which g++`  -DHDF5_BUILD_CPP_LIB:BOOL=ON  -DCMAKE_INSTALL_PREFIX:PATH=$(PWD)/hdf5/build && \
  make -j 8 && \
  make install

hgsvg/blasr/alignment/bin/blasr: hdf5/build/lib/libhdf5_cpp.a
	cd hgsvg && make  HDF5INCLUDEDIR=$(PWD)/hdf5/build/include HDF5LIBDIR=$(PWD)/hdf5/build/lib

hgsvg/blasr/pbihdfutils/bin/samtobas: hdf5/build/lib/libhdf5_cpp.a
	cd hgsvg && make

pbsamstream/pbsamstream:
	cd pbsamstream && make

samtools/samtools: local_assembly/pbgreedyphase/partitionByPhasedSNVs
	cd samtools && make CFLAGS=-I$(abspath local_assembly/pbgreedyphase/bzip2-1.0.6) LDFLAGS="-L$(abspath local_assembly/pbgreedyphase/bzip2-1.0.6) -lbz2"
