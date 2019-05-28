#!/usr/bin/env bash
echo "$1 --vcf $2 --sample $3 --region $4 --minSites $5"
res=`$1 --vcf $2 --sample $3 --region $4 --minSites $5` 
if [ $res -eq "1" ]; then
		grep -v "^@" $6 >> $7
		grep -v "^@" $6 >> $8
fi
