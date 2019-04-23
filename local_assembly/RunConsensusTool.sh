#!/usr/bin/env bash
if [  "$4" = "HGSVG_BAM" ]; then
    quiver  -j4 --minCoverage 7 --noEvidenceConsensusCall nocall --referenceFilename $1 $2 -o $3
else
    echo "Running arrow"
    arrow  -j4 --minCoverage 7 --noEvidenceConsensusCall nocall --referenceFilename $1 $2 -o $3
fi