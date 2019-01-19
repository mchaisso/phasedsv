#!/usr/bin/env bash

vcftools --gzvcf $1 --indv $2  --recode  --stdout | vcffilter -g "( GT = 1|0 | GT = 0|1 ) "  | awk '{ if ($NF != ".") print;}' | bedtools cluster -i stdin -d $3 | bedtools groupby -g 11 -c 11,1,2,2 -o count,first,first,last | awk '{ print $3"\t"$4"\t"$5"\t"$5-$4"\t"$2;}' 
