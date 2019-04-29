#!/usr/bin/env python
import sys
minLength = int(sys.argv[1])

for line in sys.stdin:
    if line[0] == '@':
        sys.stdout.write(line)
        continue
    vals = line.split()
    title = vals[0].split("/")


    srLen=len(vals[9])
    if srLen > minLength:
        sys.stdout.write(line)
            
