all: local_assembly/shiftSamPos \
	mcutils/src/samToBed \
  local_assembly/pbgreedyphase/partitionByPhasedSNVs \
  hdf5/build/lib/libhdf5_cpp.a \
  hgsvg/blasr/alignment/bin/blasr \
  hgsvg/blasr/pbihdfutils/bin/samtobas \
  pbsamstream/pbsamstream \
  samtools/samtools \
  environments/python2.7/bin/activate \
  setup_phasedsv.sh

/environments/python2.7/bin/activate:
	source ./setup_virtualenv.sh

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
	cd local_assembly/pbgreedyphase && make

hdf5/hdf5-1.8.20/README.txt:
	mkdir -p hdf5
	cd hdf5 && wget https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.20.tar.gz
	cd hdf5 && tar xvf hdf5-1.8.20.tar.gz

hdf5/build/lib/libhdf5_cpp.a: hdf5/hdf5-1.8.20/README.txt
	cd hdf5/hdf5-1.8.20 &&./configure --prefix=$(PWD)/hdf5/build --enable-cxx &&\
	make -j 8 && make install

hgsvg/blasr/alignment/bin/blasr: hdf5/build/lib/libhdf5_cpp.a
	cd hgsvg && make  HDF5INCLUDEDIR=$(PWD)/hdf5/build/include HDF5LIBDIR=$(PWD)/hdf5/build/lib

hgsvg/blasr/pbihdfutils/bin/samtobas: hdf5/build/lib/libhdf5_cpp.a
	cd hgsvg && make

pbsamstream/pbsamstream:
	cd pbsamstream && make

samtools/samtools: local_assembly/pbgreedyphase/partitionByPhasedSNVs
	cd samtools && make CFLAGS=-I$(abspath local_assembly/pbgreedyphase/bzip2-1.0.6) LDFLAGS="-L$(abspath local_assembly/pbgreedyphase/bzip2-1.0.6) -lbz2"
