#!/usr/bin/env  bash
base=dirname $0
#
#$ -t 1-4896 -tc 200 -S /bin/bash -V  -e /dev/null -o /dev/null -l mfree=2G -l h_rt=02:00:00 -pe serial 6 -wd /net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/AutozygousPartitioning/HG00733 -p -600 -q eichler-short.q -l disk_free=2G
#
$base/RunTrioTiledAssemblyOnRegions.sh `awk "NR == $SGE_TASK_ID" regions.txt ` lasm /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/local_assembly/RunTrioPartionedAssemblies.mak  HG00733.HGSVG.params.txt $SGE_TASK_ID 1

