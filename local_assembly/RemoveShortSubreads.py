#!/usr/bin/env python
import sys
minLength = int(sys.argv[1])

for line in sys.stdin:
    if line[0] == '@':
        sys.stdout.write(line)
        continue
    vals = line.split()
    title = vals[0].split("/")
    
    if len(title) == 3:
        sr = title[2].split("_")
        srLen = int(sr[1]) - int(sr[0])
        if srLen > minLength:
            sys.stdout.write(line)
            
