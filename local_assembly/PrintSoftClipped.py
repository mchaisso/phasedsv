#!/usr/bin/env python

import pysam
import argparse

ap = argparse.ArgumentParser(description="Print soft clipped reads/assemblies")
ap.add_argument("--bam", help="Input bam file", required=True)
ap.add_argument("--minClip", help="Minimum clip to report.", required=True, type=int)
ap.add_argument("--maxAlignLength", help="Do not consider a contig clipped if it is of this length or longer, NULL allows any length", type=int, default=None)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()


bamFile = pysam.AlignmentFile(args.bam, 'rb')
outFile = open(args.out, 'w')
for read in bamFile:
    clipStart = read.query_alignment_start
    clipEnd   = read.query_length - read.query_alignment_end
    alignLength = read.query_alignment_end - read.query_alignment_start
    tStart    = read.reference_start
    tEnd      = read.reference_end

    if args.maxAlignLength is not None and alignLength > args.maxAlignLength:
        continue
    if clipStart > args.minClip:
        outFile.write(read.reference_name + "\t" + str(tStart) + "\t" + str(int(tStart) + 1) + "\t" + read.query_name + "\t" + "3p" + "\t" + str(clipStart) + "\n")

    if clipEnd > args.minClip:
        outFile.write(read.reference_name + "\t" + str(tEnd) + "\t" + str(int(tEnd) + 1) + "\t" + read.query_name + "\t" + "5p" + "\t" + str(clipEnd) + "\n")

            


