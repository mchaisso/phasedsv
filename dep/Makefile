# Build PhasedSV dependencies, mirrored from SMRT-SV distribution

# Define
MINICONDA_INSTALLER = Miniconda3-4.5.11-Linux-x86_64.sh

# All
.PHONY: all
build/install_flags/all: \
	build/install_flags/env_python2_install \
	build/install_flags/env_python3_install \
	build/install_flags/env_pacbio_install \
	build/install_flags/env_tools_install \
  install_flags/dep_conda_link


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
	. bin/activate python3; \
	source $(abspath install_python3.sh); \
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

# Conda
install_flags/dep_conda_link: install_flags/dep_conda_build
	mkdir -p bin install_flags

# python3
	ln -sf ../conda/build/envs/python3/bin/python3 bin/python3
	ln -sf ../conda/build/envs/python3/bin/python3 bin/python3.6
	ln -sf ../conda/build/envs/python3/bin/ipython3 bin/ipython3
	ln -sf ../conda/build/envs/python3/bin/snakemake bin/snakemake

# python2
	ln -sf ../conda/build/envs/python2/bin/python2 bin/python2

# pacbio
	ln -sf ../conda/build/envs/pacbio/bin/python2 bin/pb_python2
	ln -sf ../conda/build/envs/pacbio/bin/blasr bin/blasr
	ln -sf ../conda/build/envs/pacbio/bin/pbindex bin/pbindex
	ln -sf ../conda/build/envs/pacbio/bin/sawriter bin/sawriter
	ln -sf ../conda/build/envs/pacbio/bin/variantCaller bin/variantCaller
	ln -sf ../conda/build/envs/pacbio/bin/quiver bin/quiver
	ln -sf ../conda/build/envs/pacbio/bin/arrow bin/arrow
	ln -sf ../conda/build/envs/pacbio/bin/bax2bam bin/bax2bam
  ln -sf ../conda/build/envs/pacbio/bin/pbmm2 bin/pbmm2

# Tools
	ln -sf ../conda/build/envs/tools/bin/canu bin/canu
	ln -sf ../conda/build/envs/tools/bin/bedtools bin/bedtools
	ln -sf ../conda/build/envs/tools/bin/groupBy bin/groupBy
  ln -sf ../conda/build/envs/tools/bin/vt bin/vt
	ln -sf ../conda/build/envs/tools/bin/vcfkeepinfo bin/vcfkeepinfo
	ln -sf ../conda/build/envs/tools/bin/vcf2tsv bin/vcf2tsv
	ln -sf ../conda/build/envs/tools/bin/vcffixup bin/vcffixup
	ln -sf ../conda/build/envs/tools/bin/tabix bin/tabix
	ln -sf ../conda/build/envs/tools/bin/bgzip bin/bgzip
	ln -sf ../conda/build/envs/tools/bin/seqtk bin/seqtk
	ln -sf ../conda/build/envs/tools/bin/samtools bin/samtools
	ln -sf ../conda/build/envs/tools/bin/bamleftalign bin/bamleftalign
	ln -sf ../conda/build/envs/tools/bin/RepeatMasker bin/RepeatMasker
	ln -sf ../conda/build/envs/tools/bin/blastp bin/blastp
	ln -sf ../conda/build/envs/tools/bin/trf bin/trf
	ln -sf ../conda/build/envs/tools/bin/bcftools bin/bcftools
	ln -sf ../conda/build/envs/tools/bin/bioawk bin/bioawk

	#mkdir -p share
	#ln -sf ../conda/build/envs/tools/share/RepeatMasker share/

	# Flag complete
	date > $@

install_flags/dep_conda_build:
	mkdir -p install_flags
	make -C conda
	date > $@


# Clean
.PHONY: clean
clean:
	make -C conda clean
	make -C smrtsvtools clean
	rm -rf bin
	rm -rf install_flags
