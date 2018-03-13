#!/usr/bin/env python

import sys
import argparse
from sets import Set
ap = argparse.ArgumentParser(description="Format fasta for falcon reading")
ap.add_argument("--fakename", help="Copy the name of the first movie into all for local assembly.", default=False, action='store_true')
ap.add_argument("--unique", help="When there are multiple reads with the same name, only print one", action='store_true', default=False)
ap.add_argument("--maxLength", help="Maximum length read to allow.", type=int, default=0)
args = ap.parse_args()

movieIndex = 1
firstMovieName = None
seen = Set()

for line in sys.stdin.readlines():
    if (line[0] == '>'):
        title = line
        origTitle = title
    else:
        if (args.fakename):
            if (firstMovieName is None):
                vals = title.rstrip().split("/")
                firstMovieName = vals[0][1:]

            title = ">{}/{}/1_{}".format(firstMovieName, movieIndex, len(line))
        movieIndex +=1

        L=50
        i = 0
        if (origTitle in seen and args.unique):
            continue
        print title
        while (i < len(line) and (args.maxLength == 0 or i < args.maxLength)):
            e = min(i+L,len(line))
            print line[i:e].strip()
            i+=L
        seen.add(origTitle)
