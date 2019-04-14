all: shiftSamPos pbgreedyphase/partitionByPhasedSNVs 

MAKE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

shiftSamPos: ShiftSamFromRegion.cpp
	g++ $< -o $@

pbgreedyphase/partitionByPhasedSNVs:
	cd pbgreedyphase && make
