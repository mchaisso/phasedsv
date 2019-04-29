#!/usr/bin/env python

import argparse

ap = argparse.ArgumentParser(description="")
ap.add_argument("--minSubreadLength", help="Minimum subread length to print", default=None, type=int)

args = ap.parse_args()

import sys
seen = {}

for line in sys.stdin:
        if line[0] == "@":
            sys.stdout.write(line)
            continue
        v = line.split()
        if v[0] in seen:
            continue
        if args.minSubreadLength is not None:
            name = v[0]
            sr = name.split("/")[2].split("_")
            subreadLen = int(sr[1]) - int(sr[0])
            if subreadLen < args.minSubreadLength:
                continue
        seen[v[0]] = True
        sys.stdout.write(">" + v[0]+ "\n")
        sys.stdout.write(v[9] +"\n")

