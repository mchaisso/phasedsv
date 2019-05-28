


BASE=$(dirname $0)
SAM=$1
READ_SOURCE=$2
ASSEMBLY=$3
BAM=$4

pbmm2 index $ASSEMBLY $ASSEMBLY.mmi
a=$(basename $SAM)
b=${a%.*}
if [[ ! -z $READ_SOURCE  && "$READ_SOURCE"=="HGSVG_BAM" ]]; then
    echo "cat $SAM | $BASE/ReformatHGSVGSam.py | samtools view -bhS - -o $b.bam"
    cat $SAM | $BASE/ReformatHGSVGSam.py | samtools view -bS - | samtools sort -o $b.bam
		pbindex $b.bam
else
    samtools view -bhS $SAM -o $b.bam
fi
pbmm2 align  $ASSEMBLY.mmi $b.bam $BAM --sort
pbindex $BAM


		

