MAKE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

all: local_assembly/shiftSamPos \
  local_assembly/pbgreedyphase/partitionByPhasedSNVs \
  install_flags/dep_conda_link \
  config.sh
SHELL := /bin/bash

install_flags/dep_conda_link:
	cd dep && make

config.sh:
	echo "export PATH=\$$PATH:"${MAKE_DIR}"/dep/bin " > $@
	echo "# Add any additional configuration here. If your system" >> $@
	echo "# uses modules, this will likely involve loading required" >> $@
	echo "# modules.">> $@

local_assembly/shiftSamPos: install_flags/dep_conda_link
	cd local_assembly && make CONDA_PREFIX=${MAKE_DIR}/dep

local_assembly/pbgreedyphase/partitionByPhasedSNVs: install_flags/dep_conda_link
	cd local_assembly && make

