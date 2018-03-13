all: shiftSamPos mcutils/src/samToBed

shiftSamPos: ShiftSamFromRegion.cpp
	g++ $< -o $@


mcutils/src/samToBed:
	cd mcutils/src && make

pbgreedyphase/partitionByPhasedSNVs:
  cd pbgreedyphase && make
