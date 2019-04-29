#!/usr/bin/env python
import sys
anyRG=None
firstNonComment=True
for line in sys.stdin:
    if len(line) >= 3 and line[0:3] == "@HD":
        line="@HD\tVN:1.MC.rc43\tpb:3.0.1\tSO:coordinate\n"
    elif len(line) >= 3 and line[0:3] == "@RG":
        line = line.strip() + ";FRAMERATEHZ=75.000000" + "\n"
        vals = line.split()
        sys.stderr.write(vals[1] + "\n")
        if vals[1][0:2]=="ID":
            if "-" in vals[1]:
                sys.stderr.write("Skipping rg " + vals[1] + "\n")
                continue
            else:
                anyRG=vals[1][4:]
    elif line[0] != "@":
        if firstNonComment == True:
            sys.stdout.write("@RG\tID:123456\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=100372700;SEQUENCINGKIT=100612400;BASECALLERVERSION=2.3.0.3.154799;FRAMERATEHZ=75.005768\tPU:m160317_220824_42274_c100962212550000001823220307011674_s1_p0\tPM:RS\n")
            firstNonComment = False
        vals=line.strip().split()

        kvps=vals[11:]
        QS=-1
        QE=-1
        qsi=-1
        qei=-1
        foundRG=False
        for i in range(0,len(kvps)):
            if kvps[i][0:5] == "sq:i:":
                qsi=i
            if kvps[i][0:5] == "qe:i:":
                qei=i
            if kvps[i][0:5] == "RG:Z:":
                kvps[i]="RG:Z:123456"

        if foundRG is False:
            kvps.append("RG:Z:123456")
        if qsi != -1 and qei != -1:
            qs=int(kvps[qsi][5:])
            qe=int(kvps[qei][5:])
            kvps[qsi] = "qs:i:0"
            kvps[qei] = "qe:i:"+ str(qe-qs)
        vals =vals[0:11] + kvps
        line = "\t".join(vals) + "\n"
    sys.stdout.write(line)
        
        
        
        
