#!/usr/bin/env python
import sys
for line in sys.stdin:
    if line[0] == '@':
        sys.stdout.write(line)
        continue
    vals = line.split()
    
