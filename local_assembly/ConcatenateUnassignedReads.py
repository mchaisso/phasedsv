#!/usr/bin/env python


import pysam
import argparse
import re
import sys
import intervaltree

ap = argparse.ArgumentParser(description="Determine how much of a region is phased")
ap.add_argument("--vcf", help="VCF file, should have trio", required=True)
ap.add_argument("--sample", help="SampleID  to consider.", required=True)
ap.add_argument("--region", help="Region to consider", required=True)
ap.add_argument("--minSites", help="If there are no more than this many sites, concatenate reads.", type=int)

args=ap.parse_args()

#
# Turn the region

vcfFile = pysam.VariantFile(args.vcf)  # auto-detect input format

regionRe = re.compile("(.*)[\.:](.*)-(.*)")

def GetBounds(a):
    if len(a) == 0:
        return 0,0
    else:
        return min(a), max(a)
    

m = regionRe.match(args.region)
g = m.groups()
chrom = g[0]
start = int(g[1])
end   = int(g[2])

auto = intervaltree.IntervalTree()
auto.addi(start, end)
    
# Fetch sites
phased = []
het    = []
hom    = []

for rec in vcfFile.fetch(chrom,start,end):
    #
    # Father is homozygous, motheris het, so it is possible to determine hertiance from mother
    #
    if rec.samples[args.sample]['PS'][0] != '.':
        phased.append(rec.start)


if len(phased) < args.minSites:
    print("1")
else:
    print("0")


