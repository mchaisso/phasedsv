#!/usr/bin/env python
import sys
import argparse
import os

ap=argparse.ArgumentParser(description="Prepare submission scripts for sge, uge, or slurm cluster management systems.")

swd=os.path.dirname(os.path.realpath(__file__))

ap.add_argument("--regions", help="Full path to regions.", required=True)
ap.add_argument("--params", help="Path to parameter file", required=True)
ap.add_argument("--runTrio", help="Prepare commands for trio assembly.", default=False, action='store_true')
ap.add_argument("--grid", help="Grid type", choices=["sge", "uge", "slurm"], required=True)
ap.add_argument("--base", help="base name for grid commands", default="psv_grid")
ap.add_argument("--conc", help="Number of concurrent jobs", default="50")
ap.add_argument("--config", help="Extra configuration parameters for job submission", default="")

args=ap.parse_args()

regionsFile = open(args.regions)
numRegions = len(regionsFile.readlines())
assemblyScript=swd+"/../AssembleRegion.sh"
if args.runTrio:
    assemblyScript=swd+"/../RunTrioTiledAssemblyOnRegions.sh"


if args.grid == "slurm":
    gridFileName=args.base+".run.sh"
    gf=open(gridFileName,'w')
    gf.write("#!/usr/bin/env bash\n")
    gf.write(assemblyScript +" `awk \"NR == $SLURM_ARRAY_TASK_ID\" " + args.regions + " ` " + args.params + "\n")
    gf.close()

    submit=args.base+".submit.sh"
    sf=open(submit,'w')
    sf.write("#!/usr/bin/env bash\n")
    sf.write("sbatch --cpus-per-task=4 --mem=4G --time=1:00:00 --array=1-{}%{} {} {}\n".format(numRegions, args.conc, gridFileName, args.config))
    sf.close()

    print "\nCreated job script '" + gridFileName + "'"
    print "Created submission script '" + submit+ "'.  You may need to modify"
    print "the submission script for custom array parameters at your site."
    

