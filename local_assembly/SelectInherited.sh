#!/usr/bin/env bash
BASE=$(dirname $0)

usage()
{
cat << EOF
usage: $0 options


OPTIONS:
   -h      Show this message
   -b      BAM-FOFN
   -v      Trio vcf
   -i      Inheritance file for sample
   -r      Region
   -m      Member (fa or mo)
   -R      Reference
   -s      Sample
EOF
}
if [ "$#" -lt 7 ]; then
		usage
		exit
fi
while getopts “hb:r:v:i:R:m:s:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         b)
             bams=$OPTARG
             ;;
         v)
             vcf=$OPTARG
             ;;
         i)
             inherited=$OPTARG
             ;;
         r)
             region=$OPTARG
             ;;
				 m)
						 member=$OPTARG
						 ;;
				 R)
						 ref=$OPTARG
						 ;;
				 s)
						 sample=$OPTARG
						 ;;
         ?)
             usage
             exit
             ;;
     esac
done


echo $region | tr ":-" "\t\t" > $member.query.bed

head -1 $bams | xargs -i samtools view -H {} > $member.reads.sam
cat $bams | xargs -i samtools view -q 30 {} $region  |  sed -e "s/qi/iq/" -e "s/qd/dq/" -e "s/qs/sq/" -e "s/qm/mq/" -e "s/td/dt/" -e "s/ts/st/" >> $member.reads.sam
$BASE/../mcutils/src/samToBed $member.reads.sam | $BASE/DetectChimeras.py > $member.filter.list


inherit=`tabix $inherited $region | bedtools intersect -a stdin -b $member.query.bed | cut -f 4`
chrom=`tabix $inherited $region | bedtools intersect -a stdin -b $member.query.bed | cut -f 5`
echo $inherit > $member.hap.txt
echo $chrom > $member.chrom.txt
nhap=`tabix $inherited $region | bedtools intersect -a stdin -b $member.query.bed | cut -f 4 | wc -l`
tabix -h $vcf $region > $member.vcf
if [ $nhap == 1 ]; then
		echo "$BASE/pbgreedyphase/partitionByPhasedSNVs --vcf $member.vcf --sam $member.reads.sam --rgn $region --pad 100000 --h1 $member.h0.sam --h2 $member.h1.sam --ref $ref --minGenotyped 1 --summary summary.txt --sample $sample"
		grep -v -f filter.list $member.reads.sam | $BASE/RemoveShortSubreads.py 1000 | $BASE/pbgreedyphase/partitionByPhasedSNVs --vcf $member.vcf --sam /dev/stdin --rgn $region --pad 100000 --h1 $member.h0.sam --h2 $member.h1.sam --ref $ref --minGenotyped 1 --summary summary.txt --sample $sample --unassigned /dev/null
		
		if [ $inherit == 0 ]; then
				cp $member.h0.sam $member.inherited.sam
				
		fi
		if [ $inherit == 1 ]; then
				cp $member.h1.sam $member.inherited.sam
		fi
fi


touch $member.inherited.sam

		

