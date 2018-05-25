#!/usr/bin/env bash

usage()
{
cat << EOF
		RunTrioTiledAssemblyOnRegions.sh regions paramfile [ -j jobname ] [-i idx ] [ -d delay]
		-j job  Will submit jobs under this name to sge.
		-i idx  The index of this job, 1... #jobs. When less than delay, sleep
		               for a few seconds to avoid slamming a shared filesystem with 
		               too many jobs that read from bams at once. It is assumed that 
		               after a certain number of jobs the sysem will be in steady state.
		               Typically, with 200 jobs submitted, steady state is ~1 new job every
		               3s.
		-d del  Delay job starts for indices less than idx
EOF
exit 1

}
if [ -z $TMPDIR  ]; then
		echo "ERROR. TMPDIR must be set, and should be a full path from /."
		echo "The assembly pipeline will execut as a subdirectory of  \$TMPDIR."
		exit 1
fi

BASE="$( cd "$(dirname "$0")" ; pwd -P )"


if [ $# -lt 2 ]; then
		usage
		exit
fi
REGIONS=$1
shift
PARAMFILE=$1
shift
JOBNAME="loc_asm"
DELAY=1
INDEX=1
while getopts “hj:i:d:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         j)
             JOBNAME=$OPTARG
             ;;
         i)
             INDEX=$OPTARG
             ;;
         d)
             DELAY=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done


INDEX=1
DODELAY=1
JOBNAME=loc_asm

ASSEMBLER=$BASE/RunTrioPartionedAssemblies.mak




if [ $INDEX -lt $DODELAY ]; then
		v=$(($INDEX*5))
		sleep $v"s"
fi

DIR=$TMPDIR/$USER/lasm/$$/
echo "starting" $DIR



PARAMS=`cat $PARAMFILE | tr "\n" " "`
source $PARAMFILE
source $BASE/../setup_phasedsv.sh
echo $LD_LIBRARY_PATH

mkdir -p $DEST/samfiles
mkdir -p $DEST/assemblies
mkdir -p $DEST/records

mkdir -p $DIR

echo $REGIONS | tr ";" "\n" | sed '/^\s$/d'>  $DIR/Regions.txt


cat $DIR/Regions.txt | tr "\.-" "\t\t" > $DIR/Regions.bed
bedtools sort -i $DIR/Regions.bed | bedtools merge | awk '{ print $1":"$2"-"$3;}' > $DIR/Regions.rgn


for r in `cat $DIR/Regions.rgn`; do
		n=`echo $r | tr ":" "."`
		echo $BAMS
		samtools merge -u -f -b $BAMS -R $r $DIR/reads.region.ch.$n.bam
		samtools index $DIR/reads.region.ch.$n.bam
		echo $MOBAMS
		samtools merge -u -f -b $MOBAMS -R $r $DIR/reads.region.mo.$n.bam
		samtools index $DIR/reads.region.mo.$n.bam
		echo $FABAMS
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

echo "Done configuring fofn"

#
# Set up param file for local files.
#
grep -v "BAMS=" $PARAMFILE > $DIR/params.txt
grep "DEST" $PARAMFILE >> $DIR/params.txt
echo "BAMS=$DIR/reads.region.ch.bam.fofn" >> $DIR/params.txt
echo "MOBAMS=$DIR/reads.region.mo.bam.fofn" >> $DIR/params.txt
echo "FABAMS=$DIR/reads.region.fa.bam.fofn" >> $DIR/params.txt
echo "cat $DIR/Regions.txt | xargs -P 3 -I {} $BASE/RunTiledAssembly.sh {} $JOBNAME $BASE/RunTrioPartionedAssemblies.mak  $DIR/params.txt $DIR "

#
# Local (to node) assemble cluster.
#
cat $DIR/Regions.txt | xargs -P 3 -I {} $BASE/RunTiledAssembly.sh {} $DIR/params.txt -j $JOBNAME -a $BASE/RunTrioPartionedAssemblies.mak -d $DIR
