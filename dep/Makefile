# Build PhasedSV dependencies, mirrored from SMRT-SV distribution

# Define
MINICONDA_INSTALLER = Miniconda3-4.5.11-Linux-x86_64.sh
MAKE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SHELL=/bin/bash

# All
.PHONY: all
build/install_flags/all: \
	build/install_flags/env_python2_install \
	build/install_flags/env_python3_install \
	build/install_flags/env_pacbio_install \
	build/install_flags/env_tools_install \
        build/install_flags/dep_mcutils \
        install_flags/dep_conda_link

#        build/install_flags/shift_sam \
# python: Refuse unambiguous "python" interpreter (require "python2" or "python3")
bin/python: python_dummy/python_dummy.sh
	mkdir -p bin
	cp -f $< $@
	chmod 555 $@

# Deploy and setup directories for deploying environments
build/bin/conda: build/${MINICONDA_INSTALLER}
	cd build; bash ${MINICONDA_INSTALLER} -f -b -p $(abspath build)
	mkdir -p build/install_flags
	touch $@

# Download
build/${MINICONDA_INSTALLER}:
	mkdir -p build
	wget https://repo.continuum.io/miniconda/${MINICONDA_INSTALLER} -O $@

build/install_flags/env_python2_install: build/install_flags/env_python2_deploy
	@set -e; \
	cd build; \
	. bin/activate python2; \
        echo "source $(abspath install_python2.sh)"; \
	source $(abspath install_python2.sh); \
	conda list > PACKAGES_PYTHON2; \
	. bin/deactivate
	date > $@

# Deploy environment
build/install_flags/env_python2_deploy: build/bin/conda
	@set -e; \
	cd build; \
	bin/conda create --yes --name python2 python=2.7.15
	date > $@

# Install packages
build/install_flags/env_python3_install: build/install_flags/env_python3_deploy
	@set -e; \
	cd build; \
	. bin/activate python3 ; \
        echo "MD" $(MAKE_DIR)install_python3.sh; \
	source $(MAKE_DIR)install_python3.sh; \
	conda list > PACKAGES_PYTHON3; \
	. bin/deactivate
	date > $@

# Deploy environment
build/install_flags/env_python3_deploy: build/bin/conda
	@set -e; \
	cd build; \
	bin/conda create --yes --name python3 python=3.6.2
	date > $@

# Install packages
build/install_flags/env_tools_install: build/install_flags/env_tools_deploy
	@set -e; \
	cd build; \
	. bin/activate tools; \
	source $(abspath install_tools.sh); \
	conda list > PACKAGES_TOOLS; \
	. bin/deactivate
	date > $@


# Deploy environment
build/install_flags/env_tools_deploy: build/bin/conda
	@set -e; \
	cd build; \
	bin/conda create --yes --name tools python=3.6.2
	date > $@

### Environment: pacbio - PacBio conda toolkit ###

# Install packages
build/install_flags/env_pacbio_install: build/install_flags/env_pacbio_deploy
	@set -e; \
	cd build; \
	. bin/activate pacbio; \
	source $(abspath install_pacbio.sh); \
	conda list > PACKAGES_PACBIO; \
	. bin/deactivate
	date > $@


# Deploy environment
build/install_flags/env_pacbio_deploy: build/bin/conda
	@set -e; \
	cd build; \
	bin/conda create --yes --name pacbio python=2.7.15
	date > $@

build/install_flags/dep_mcutils:
	@set -e ; \
	cd mcutils/src; \
	make
	date > $@

build/install_flags/dep_blasrmc:
	@set -e ; \
        cd blasr ; \
        make
	date > $@

build/install_flags/shift_sam:
	@set -e ; \
        cd ../local_assembly ; \
        make
	date > $@

# Conda
install_flags/dep_conda_link: build/install_flags/env_pacbio_deploy   build/install_flags/env_pacbio_install \
  build/install_flags/env_tools_deploy \
  build/install_flags/env_python3_deploy \
  build/install_flags/env_python3_install \
  build/install_flags/dep_blasrmc \
  build/install_flags/dep_mcutils

	mkdir -p bin install_flags
	mkdir -p lib
	mkdir -p include
	ln -sf $(MAKE_DIR)build/envs/pacbio/include/htslib $(MAKE_DIR)/include/htslib
	ln -sf $(MAKE_DIR)build/envs/pacbio/lib/libhts.a $(MAKE_DIR)/lib/libhts.a


# mcutils
	ln -sf $(MAKE_DIR)mcutils/src/samToBed $(MAKE_DIR)/bin/samToBed
	ln -sf $(MAKE_DIR)mcutils/src/samLiftover $(MAKE_DIR)/bin/samLiftover

# blasr
	ln -sf $(MAKE_DIR)blasr/alignment/bin/blasrmc $(MAKE_DIR)/bin/blasrmc

# python3
	ln -sf $(MAKE_DIR)build/envs/python3/bin/python3 $(MAKE_DIR)/bin/python3
	ln -sf $(MAKE_DIR)build/envs/python3/bin/python3 $(MAKE_DIR)/bin/python3.6
	ln -sf $(MAKE_DIR)build/envs/python3/bin/ipython3 $(MAKE_DIR)/bin/ipython3
	ln -sf $(MAKE_DIR)build/envs/python3/bin/snakemake $(MAKE_DIR)/bin/snakemake

# python2
	ln -sf $(MAKE_DIR)build/envs/python2/bin/python2 $(MAKE_DIR)/bin/python2
	ln -sf $(MAKE_DIR)build/envs/python2/bin/python2 $(MAKE_DIR)/bin/python

# pacbio
	ln -sf $(MAKE_DIR)build/envs/pacbio/bin/python2 $(MAKE_DIR)/bin/pb_python2
	ln -sf $(MAKE_DIR)build/envs/pacbio/bin/blasr $(MAKE_DIR)/bin/blasr
	ln -sf $(MAKE_DIR)build/envs/pacbio/bin/pbindex $(MAKE_DIR)/bin/pbindex
	ln -sf $(MAKE_DIR)build/envs/pacbio/bin/sawriter $(MAKE_DIR)/bin/sawriter
	ln -sf $(MAKE_DIR)build/envs/pacbio/bin/variantCaller $(MAKE_DIR)/bin/variantCaller
	ln -sf $(MAKE_DIR)build/envs/pacbio/bin/quiver $(MAKE_DIR)/bin/quiver
	ln -sf $(MAKE_DIR)build/envs/pacbio/bin/arrow $(MAKE_DIR)/bin/arrow
	ln -sf $(MAKE_DIR)build/envs/pacbio/bin/bax2bam $(MAKE_DIR)/bin/bax2bam
	ln -sf $(MAKE_DIR)build/envs/pacbio/bin/pbmm2 $(MAKE_DIR)/bin/pbmm2

# Tools
	ln -sf $(MAKE_DIR)build/envs/tools/bin/canu $(MAKE_DIR)/bin/canu
	ln -sf $(MAKE_DIR)build/envs/tools/bin/bedtools $(MAKE_DIR)/bin/bedtools
	ln -sf $(MAKE_DIR)build/envs/tools/bin/vt $(MAKE_DIR)/bin/vt
	ln -sf $(MAKE_DIR)build/envs/tools/bin/tabix $(MAKE_DIR)/bin/tabix
	ln -sf $(MAKE_DIR)build/envs/tools/bin/bgzip $(MAKE_DIR)/bin/bgzip
	ln -sf $(MAKE_DIR)build/envs/tools/bin/samtools $(MAKE_DIR)/bin/samtools
	ln -sf $(MAKE_DIR)build/envs/tools/bin/RepeatMasker $(MAKE_DIR)/bin/RepeatMasker
	ln -sf $(MAKE_DIR)build/envs/tools/bin/trf $(MAKE_DIR)/bin/trf
	ln -sf $(MAKE_DIR)build/envs/tools/bin/bioawk $(MAKE_DIR)/bin/bioawk
	ln -sf $(MAKE_DIR)build/envs/tools/bin/bedToBigBed $(MAKE_DIR)/bin/bedToBigBed

	# Flag complete
	date > $@


# Clean
.PHONY: clean
clean:
	rm -rf bin
	rm -rf install_flags
