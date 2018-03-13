#!/usr/bin/env python

import sys
import collections


import argparse
ap = argparse.ArgumentParser(description="Filter reads from sam")
ap.add_argument("--input", help="Input file", default="/dev/stdin")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

reads = collections.defaultdict(list)


inFile = open(args.input)
outFile = open(args.out, 'w')
for line in inFile:
    v = line.split()
    reads[v[3]].append((v[0], int(v[1]), int(v[2]), int(v[4]), int(v[5]), int(v[6])))


def Overlaps(a,b, idx=1):
    if (a[0] == b[0] and 
        ((a[idx] <= b[idx] and a[1+idx]  > b[idx]) or
         (a[idx] < b[1+idx] and a[1+idx] >= b[1+idx]))):
        return True

def Duplicate(a,b):
    return (a==b)

for readName, readAlns in reads.iteritems():
    if (len(readAlns) > 1):
        for i in range(0,len(readAlns)-1):
            for j in range(i+1,len(readAlns)):
                if (Overlaps(readAlns[i], readAlns[j]) and
                    Overlaps(readAlns[i], readAlns[j], 4) and
                    readAlns[i][3] != readAlns[j][3]):
                    outFile.write( readName + "\n")
                if (Duplicate(readAlns[i], readAlns[j])):
                    outFile.write(readName + "\n")
                    
                
            
        
    
