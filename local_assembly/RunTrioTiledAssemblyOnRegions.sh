#!/usr/bin/env bash
REGIONS=$1
JOBNAME=$2
ASSEMBLER=$3
PARAMFILE=$4
INDEX=$5
DODELAY=$6


if [ $INDEX -lt $DODELAY ]; then
		v=$(($INDEX*5))
		sleep $v"s"
fi
echo "starting"
DIR=/$TMPDIR/$USER/lasm/$$/
mkdir -p samfiles
mkdir -p assemblies
mkdir -p records



PARAMS=`cat $PARAMFILE | tr "\n" " "`
source $PARAMFILE


source /etc/profile.d/modules.sh
module load mpc/0.8.2; module load mpfr/3.1.0; module load gmp/5.0.2; module load gcc/latest
module load python/2.7.3 ;
module load htslib/1.3.2
module load samtools/latest
module load bedtools/latest

mkdir -p $DIR

echo $REGIONS | tr ";" "\n" | sed '/^\s$/d'>  $DIR/Regions.txt


cat $DIR/Regions.txt | tr "\.-" "\t\t" > $DIR/Regions.bed
bedtools sort -i $DIR/Regions.bed | bedtools merge | awk '{ print $1":"$2"-"$3;}' > $DIR/Regions.rgn


for r in `cat $DIR/Regions.rgn`; do
		n=`echo $r | tr ":" "."`
		samtools merge -u -f -b $CHBAMS -R $r $DIR/reads.region.ch.$n.bam
		samtools index $DIR/reads.region.ch.$n.bam

		samtools merge -u -f -b $MOBAMS -R $r $DIR/reads.region.mo.$n.bam
		samtools index $DIR/reads.region.mo.$n.bam

		samtools merge -u -f -b $FABAMS -R $r $DIR/reads.region.fa.$n.bam
		samtools index $DIR/reads.region.fa.$n.bam

done
echo "DONE collecting data from NFS"
samtools merge -u -f $DIR/reads.region.ch.bam $DIR/reads.region.ch.*.bam
samtools merge -u -f $DIR/reads.region.mo.bam $DIR/reads.region.mo.*.bam
samtools merge -u -f $DIR/reads.region.fa.bam $DIR/reads.region.fa.*.bam

samtools index $DIR/reads.region.ch.bam
samtools index $DIR/reads.region.mo.bam
samtools index $DIR/reads.region.fa.bam

ls $DIR/reads.region.ch.bam > $DIR/reads.region.ch.bam.fofn
ls $DIR/reads.region.mo.bam > $DIR/reads.region.mo.bam.fofn
ls $DIR/reads.region.fa.bam > $DIR/reads.region.fa.bam.fofn


grep -v "CHBAMS=" $PARAMFILE > $DIR/params.txt
echo "CHBAMS=$DIR/reads.region.ch.bam.fofn" >> $DIR/params.txt
echo "MOBAMS=$DIR/reads.region.mo.bam.fofn" >> $DIR/params.txt
echo "FABAMS=$DIR/reads.region.fa.bam.fofn" >> $DIR/params.txt

echo "cat $DIR/Regions.txt | xargs -P 3 -I {} /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/local_assembly/RunTiledAssembly.sh {} $JOBNAME /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/local_assembly/RunTrioPartionedAssemblies.mak  $DIR/params.txt $DIR "

cat $DIR/Regions.txt | xargs -P 3 -I {} /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/local_assembly/RunTiledAssembly.sh {} $JOBNAME /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/local_assembly/RunTrioPartionedAssemblies.mak $DIR/params.txt $DIR
