#!/usr/bin/env python
import sys
for line in sys.stdin:
    if len(line) >= 3 and line[0:3] == "@HD":
        line="@HD\tVN:1.MC.rc43\tpb:3.0.1\tSO:coordinate\n"
    elif len(line) >= 3 and line[0:3] == "@RG":
        line = line.strip() + ";FRAMERATEHZ=75.000000" + "\n"
        vals = line.split()
        if vals[1][0:2]=="ID":
            if "-" in vals[1]:
                sys.stderr.write("Skipping rg " + vals[1] + "\n")
                continue
    elif line[0] != "@":
        vals=line.strip().split()

        kvps=vals[11:]
        QS=-1
        QE=-1
        qsi=-1
        qei=-1
        for i in range(0,len(kvps)):
            if kvps[i][0:5] == "sq:i:":
                qsi=i
            if kvps[i][0:5] == "qe:i:":
                qei=i
            if kvps[i][0:5] == "RG:Z:" and "-" in kvps[i]:
                kvps[i] = kvps[i].split("-")[0]
        if qsi != -1 and qei != -1:
            qs=int(kvps[qsi][5:])
            qe=int(kvps[qei][5:])
            kvps[qsi] = "qs:i:0"
            kvps[qei] = "qe:i:"+ str(qe-qs)
        vals =vals[0:11] + kvps
        line = "\t".join(vals) + "\n"
    sys.stdout.write(line)
        
        
        
        
