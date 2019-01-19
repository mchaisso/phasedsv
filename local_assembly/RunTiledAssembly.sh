#!/usr/bin/env bash
BASE="$( cd "$(dirname "$0")" ; pwd -P )"

usage()
{
cat << EOF
		RunTrioTiledAssemblyOnRegions.sh regions paramfile 
		-j job  Will submit jobs under this name to sge.
		-a asm  The full path to the assembler makefile to use.
    -d dir  Run in this directory
EOF
exit 1

}

ORIG_REGION=$1
REGION=$1
shift
PARAMFILE=$1
shift
DIR=$TMPDIR/"asm"
DEST="asm"
ASSEMBLER=$BASE/RunPartitionedAssembly.mak

info=`echo $REGION | tr "/" "\t" | cut -f 2-`
REGION=`echo $REGION | tr "/" "\t" | cut -f 1`

AUTO="HET"
echo $ORIG_REGION | grep -q "HOM"

if [ $? -eq 0 ]; then
		AUTO="HOM"
fi

while getopts “hj:a:d:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         j)
             JOBNAME=$OPTARG
             ;;
				 a)
						 ASSEMBLER=$OPTARG
						 ;;
				 d)
						 DIR=$OPTARG
						 ;;
         ?)
             usage
             exit
             ;;
     esac
done


TARGETREGION=`echo $REGION | tr '.' ':'`
echo "target region " $TARGETREGION


PARAMS=`grep -v "^#" $PARAMFILE | tr "\n" " "`
PARAMS="$PARAMS IS_AUTO=$AUTO"

source $PARAMFILE

mkdir -p $DEST/samfiles
mkdir -p $DEST/assemblies
mkdir -p $DEST/records

#
#  Setup the environment
#

echo "setting up"
source $BASE/../setup_phasedsv.sh
if [ -e $DEST/assemblies/$REGION.fasta ]; then
	 echo $DEST/assemblies/$REGION.fasta
	 echo "exit early"
exit 0
fi

mkdir -p $DIR/$REGION;
p=$PWD
pushd $DIR/$REGION;
echo "RUNNING IN " $DIR/$REGION


echo "make -f $ASSEMBLER REGION=$TARGETREGION $PARAMS"
make -f $ASSEMBLER REGION=$TARGETREGION $PARAMS || true;


popd
mv -f $DIR/$REGION/assembly.consensus.fasta ./$DEST/assemblies/$REGION.fasta
mv -f $DIR/$REGION/assembly.consensus.fasta.sam ./$DEST/samfiles/$REGION.sam  
#mv -f $DIR/$REGION/summary.txt $DEST/records/$REGION.txt
rm -rf $DIR/$REGION

exit 0
