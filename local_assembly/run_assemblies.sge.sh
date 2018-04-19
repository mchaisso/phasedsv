#!/usr/bin/env  bash
#
#$ -t 1-148923 -tc 400 -S /bin/bash -V  -e /dev/null -o /dev/null -l mfree=3G -l h_rt=02:00:00 -pe serial 6 -p -600 
#

base=dirname $0

$base/RunTrioTiledAssemblyOnRegions.sh `awk "NR == $SGE_TASK_ID" Windows.60kb-span.20kbp-stride.txt ` assembly.params 

