#!/usr/bin/env bash

TARGETREGION=`echo $1 | tr '.' ':'`
echo $TARGETREGION

DIR=$5
mkdir -p samfiles
mkdir -p assemblies
mkdir -p records

ASSEMBLER=$3
PARAMFILE=$4

PARAMS=`cat $PARAMFILE | tr "\n" " "`
source $PARAMFILE

mkdir -p $DIR/$1;
cd $DIR/$1;
echo "RUNNING IN " $DIR/$1
#
#  Setup the environment

#
echo "setting up"
source /etc/profile.d/modules.sh
module load mpc/0.8.2; module load mpfr/3.1.0; module load gmp/5.0.2; module load gcc/latest
module load python/2.7.3 ;

if [ -e $DEST/assemblies/$1.fasta ]; then
	 echo $DEST/assemblies/$1.fasta
	 echo "exit early"
exit 0
fi

echo "running makefile"
echo make -f $ASSEMBLER REGION=$TARGETREGION $PARAMS
make -f $ASSEMBLER REGION=$TARGETREGION $PARAMS || true;
mkdir -p $DEST/assemblies
mkdir -p $DEST/samfiles
mkdir -p $DEST/results
echo "MOVING $DIR/$1/assembly.consensus.fasta $DEST/assemblies/$1.fasta"
mv -f $DIR/$1/assembly.consensus.fasta $DEST/assemblies/$1.fasta
mv -f $DIR/$1/assembly.consensus.fasta.sam $DEST/samfiles/$1.sam  
mv -f $DIR/$1/records.txt $DEST/records/$1.txt


exit 0
