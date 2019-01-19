#!/usr/bin/env python
import sys
import argparse
import os
import math
import subprocess

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
assemblyScript=swd+"/../RunTiledAssembly.sh"
if args.runTrio:
    assemblyScript=swd+"/../RunTrioTiledAssemblyOnRegions.sh"

def GetMaxTasksSLURM():
    cmd="scontrol --details show config"
    resStr=subprocess.check_output(cmd.split())

    res = str(resStr).split("\n")
    for line in res:
        vals = line.strip().split("=")
        if vals[0].strip() == "MaxArraySize":
            return int(vals[1].strip())
    return 1000

def CountLines(fn):
    f=open(fn)
    return len(f.readlines())

if args.grid == "slurm":
    maxTasks = GetMaxTasksSLURM()
    gridFileName=args.base+".run.sh"
    gf=open(gridFileName,'w')
    nLines=CountLines(args.regions)
    gf.write("#!/usr/bin/env bash\n")
    
    gf.write("for i in {0.."+str(nLines-1) +"}; do \n")
    gf.write(" let rem=\" $i % $1\";\n")
    gf.write("  if [ $rem == \"$SLURM_ARRAY_TASK_ID\" ]; then \n")
    gf.write(assemblyScript +"` awk \"NR == $i\" " + args.regions + "` " + args.params + "\n")
    gf.write("  fi\n")
    gf.write("done")
    gf.close()

    submit=args.base+".submit.sh"
    sf=open(submit,'w')
    sf.write("#!/usr/bin/env bash\n")
#    sf.write("cat {regions} | xargs -I@ -P {jobs} srun --cups-per-task=4 --mem=4G --time=1:00:00  ".format(regions=args.regions, jobs=args.conc) + args.base+".run.sh @\n")
    sf.write("#SBATCH --cpus-per-task=4\n#SBATCH --mem=4G\n#SBATCH --time=1:00:00\n#SBATCH --exclusive\n")
    sf.write("#SBATCH --array=0-{}%{}\n".format(maxTasks-1, args.conc))
    if (args.config is not None):
        sf.write("#SBATCH " + args.config.replace("\\","") + "\n")    
    sf.write(os.getcwd() + "/" + gridFileName + " " + str(maxTasks)+ "\n")
    sf.close()

    print("\nCreated job script '" + gridFileName + "'")
    print("Created submission script '" + submit+ "'.  You may need to modify")
    print("the submission script for custom array parameters at your site.")




def GetMaxTasks():
    cmd="qconf -sconf"
    res = subprocess.check_output(cmd.split()).split("\n")
    for line in res:
        vals = line.split()
        if vals[0] == "max_aj_tasks":
            return int(vals[1])
    # Expected default value
    return 75000
    
    

    

if args.grid == "sge":
    maxJobsPerFile=GetMaxTasks()

    nFiles = int(math.ceil(float(numRegions)/maxJobsPerFile))
    startJob=1
    allJobs=[]
    for i in range(0,nFiles):
        gridFileName=args.base + "." + str(i) + ".run.sh"
        gf=open(gridFileName,'w')
        gf.write("#!/usr/bin/env bash\n")
        endJob=min(numRegions, startJob+maxJobsPerFile-1)
        gf.write("#$ -t {}-{} -tc {} -S /bin/bash -V  -e /dev/null -o /dev/null -l mfree=2G -l h_rt=01:00:00 -pe serial 4\n".format(startJob, endJob, args.conc))
        gf.write(assemblyScript +" `awk \"NR == $SGE_TASK_ID\" " + args.regions + " ` " + args.params + "\n")
        startJob+=maxJobsPerFile
        gf.close()
        allJobs.append(gridFileName)
    
    submit=args.base+".submit.sh"
    sf=open(submit,'w')
    sf.write("#!/usr/bin/env bash\n")
    for i in range(0,len(allJobs)):
        if i > 0:
            hold = " -hold_jid $h{} ".format(i-1)
        else:
            hold = ""
        sf.write("output=`qsub -cwd {} {}`\n".format(hold, allJobs[i]))
        sf.write("h{}=`echo $output | tr \" \" \"\\n\"  | grep -A 1 job | tail -1 | tr \".\" \"\\t\" | cut -f 1`\n".format(i))
    sf.close()
        
    print("\nCreated job script(s) '" + ", ".join(allJobs) + "'")
    print("Created submission script '" + submit+ "'.  You may need to modify")
    print("the submission script for custom array parameters at your site.")
