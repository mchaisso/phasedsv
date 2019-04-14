all: local_assembly/shiftSamPos \
	mcutils/src/samToBed \
  local_assembly/pbgreedyphase/partitionByPhasedSNVs \
  install_flags/dep_conda_link
  config.sh
SHELL := /bin/bash

install_flags/dep_conda_link:
	cd dep && make

config.sh: setup_phasedsv.sh
	echo "source $(PWD)/setup_phasedsv.sh" > $@
	echo "# Add any additional configuration here. If your system" >> $@
	echo "# uses modules, this will likely involve loading required" >> $@
	echo "# modules.">> $@

local_assembly/shiftSamPos: install_flags/dep_conda_link
	cd local_assembly && make CONDA_PREFIX=$(PWD)/dep

local_assembly/pbgreedyphase/partitionByPhasedSNVs: install_flags/dep_conda_link
	cd local_assembly && make

mcutils/src/samToBed:
	cd mcutils/src && make -j 8

