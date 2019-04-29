#!/usr/bin/env python
import sys

infile = open(sys.argv[1])
outFile = open(sys.argv[2],'w')
minNReads=int(sys.argv[3])
minCov=float(sys.argv[4])
writeLine = False
for line in infile:
    if len(line) > 0 and line[0] == ">":
        vals=line.split(" ")
        writeLine=True
        for kvp in vals:
            kv=kvp.split("=")
            if kv[0] == "reads":
                nReads=int(kv[1])
                if nReads < minNReads:
                    writeLine=False
            if kv[0] == "covStat":
                cov=float(kv[1])
                if cov < minCov:
                    writeLine=False
    if writeLine:
        outFile.write(line)
