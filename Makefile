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
  environments/python2.7/lib/python2.7/site-packages/h5py-2.8.0.post0-py2.7-linux-x86_64.egg \
  setup_phasedsv.sh \
  config.sh
SHELL := /bin/bash

environments/python2.7/bin/activate:
	./setup_virtualenv.sh

environments/python2.7/lib/python2.7/site-packages/h5py-2.8.0.post0-py2.7-linux-x86_64.egg: environments/python2.7/bin/activate hdf5/build/lib/libhdf5_cpp.a
	source ./environments/python2.7/bin/activate && \
   cd h5py && \
   python setup.py configure --hdf5=$(PWD)/hdf5/build && \
   python setup.py build && \
   python setup.py install && \
   python setup.py install_egg_info 



quiver/lib/python2.7/site-packages/pbcommand-1.0.0-py2.7.egg:
	mkdir -p quiver/lib/python2.7/site-packages/
	cd pbcommand && python setup.py build && \
    export PYTHONPATH=$$PYTHONPATH:$(PWD)/quiver/lib/python2.7/site-packages/ && \
    python setup.py install  --prefix=$(PWD)/quiver/


setup_phasedsv.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "source "$(PWD)"/environments/python2.7/bin/activate" >> $@
	echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:"$(PWD)"/quiver/lib/">> $@
	echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:"$(PWD)"/quiver/lib64/">> $@
	echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:"$(PWD)"/hdf5/build/lib/">> $@
	echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:"$(PWD)"/local_assembly/pbgreedyphase/lzma/build/lib">> $@
	echo "export PYTHONPATH=\$$PYTHONPATH:"$(PWD)"/quiver/lib/python2.7/site-packages/">> $@
	echo "export PATH=\$$PATH:"$(PWD)"/quiver/bin/" >> $@
	echo "export PATH=\$$PATH:"$(PWD)"/bin/" >> $@	
	echo "#" >> $@
	echo "# Add custom configuration here." >> $@
	echo "#" >> $@

config.sh: setup_phasedsv.sh
	echo "source $(PWD)/setup_phasedsv.sh" > $@
	echo "# Add any additional configuration here. If your system" >> $@
	echo "# uses modules, this will likely involve loading required" >> $@
	echo "# modules.">> $@

bin/vt:
	mkdir -p bin
	cd vt && make -j 8; cp vt ../bin

local_assembly/shiftSamPos: hdf5/build/lib/libhdf5_cpp.a
	cd local_assembly && make

mcutils/src/samToBed:
	cd mcutils/src && make -j 8

local_assembly/pbgreedyphase/partitionByPhasedSNVs: hdf5/build/lib/libhdf5_cpp.a
	cd local_assembly && make

local_assembly/blasr/alignment/bin/blasr: hdf5/build/lib/libhdf5_cpp.a
	cd local_assembly && make

hdf5-1.8.14:
	wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz && tar xvf hdf5-1.8.14.tar.gz


#        -DBUILD_SHARED_LIBS:BOOL=ON -std=c++11

hdf5/build/lib/libhdf5_cpp.a: hdf5-1.8.14
	rm -rf $(PWD)/hdf5-1.8.14/cmake_build
	export CXXFLAGS="-std=c++11"
	cd hdf5-1.8.14/ && \
  mkdir cmake_build && \
  cd cmake_build && \
  CXXFLAGS="-std=c++98 -D_GLIBCXX_USE_CXX11_ABI=0" && cmake .. -DCMAKE_C_COMPILER=`which gcc` \
        -DCMAKE_CPP_COMPILER=`which g++` \
        -DHDF5_BUILD_CPP_LIB:BOOL=ON \
        -DHDF5_BUILD_HL_LIB:BOOL=ON \
        -DCMAKE_INSTALL_PREFIX:PATH=$(PWD)/hdf5/build && \
  make -j 8 VERBOSE=1 && \
  make install

hgsvg/blasr/alignment/bin/blasr: hdf5/build/lib/libhdf5_cpp.a
	cd hgsvg && make  HDF5INCLUDEDIR=$(PWD)/hdf5/build/include HDF5LIBDIR=$(PWD)/hdf5/build/lib -j 8

hgsvg/blasr/pbihdfutils/bin/samtobas: hdf5/build/lib/libhdf5_cpp.a
	cd hgsvg && make -j 8 HDF5INCLUDEDIR=$(PWD)/hdf5/build/include HDF5LIBDIR=$(PWD)/hdf5/build/lib 

pbsamstream/pbsamstream:
	cd pbsamstream && make

samtools/samtools: local_assembly/pbgreedyphase/partitionByPhasedSNVs
	cd samtools && make CFLAGS=-I$(abspath local_assembly/pbgreedyphase/bzip2-1.0.6) LDFLAGS="-L$(abspath local_assembly/pbgreedyphase/bzip2-1.0.6) -lbz2"
	mkdir -p bin
	cp samtools/samtools bin/

