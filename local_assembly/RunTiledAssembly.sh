#!/usr/bin/env bash
BASE=$(readlink -f $(dirname $0))
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

TARGETREGION=`echo $1 | tr '.' ':'`
echo $TARGETREGION


REGION=$1
shift
PARAMFILE=$1
shift
DIR="asm"
ASSEMBLER=$BASE/RunPartitionedAssembly.mak
echo "running RunTiledassembly.sh"
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


PARAMS=`grep -v "^#" $PARAMFILE | tr "\n" " "`
source $PARAMFILE

mkdir -p $DEST/samfiles
mkdir -p $DEST/assemblies
mkdir -p $DEST/records

mkdir -p $DIR/$REGION;
p=$PWD
pushd $DIR/$REGION;
echo "RUNNING IN " $DIR/$REGION
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

echo "running makefile"
echo make -f $ASSEMBLER REGION=$TARGETREGION $PARAMS
make -f $ASSEMBLER REGION=$TARGETREGION $PARAMS || true;

echo "MOVING $DIR/$REGION/assembly.consensus.fasta $DEST/assemblies/$REGION.fasta"
popd
mv -f $DIR/$REGION/assembly.consensus.fasta ./$DEST/assemblies/$REGION.fasta
mv -f $DIR/$REGION/assembly.consensus.fasta.sam ./$DEST/samfiles/$REGION.sam  
mv -f $DIR/$REGION/summary.txt $DEST/records/$REGION.txt


exit 0
