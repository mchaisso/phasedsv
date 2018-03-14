all: local_assembly/shiftSamPos \
  local_assembly/pbgreedyphase/partitionByPhasedSNVs \
  lib/libhdf5_cpp.a \
  hgsvg/blasr/alignment/bin/blasr \
  pbsamstream/pbsamstream


local_assembly/shiftSamPos:
	cd local_assembly && make


local_assembly/pbgreedyphase/partitionByPhasedSNVs:
	cd local_assembly/pbgreedyphase && make

hdf5/hdf5-1.8.20/README.txt:
	mkdir -p hdf5
	cd hdf5 && wget https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.20.tar.gz
	cd hdf5 && tar xvf hdf5-1.8.20.tar.gz

hdf5/build/lib/libhdf5_cpp.a: hdf5/hdf5-1.8.20/README.txt
	cd hdf5/hdf5-1.8.20 &&./configure --prefix=$(PWD)/hdf5/build --enable-cxx &&\
    make -j 8 && make install

hgsvg/blasr/alignment/bin/blasr:
	cd hgsvg && make

pbsamstream/pbsamstream:
	cd pbsamstream && make







