#!/usr/bin/env python
import sys
import argparse
import os
import math
import subprocess

ap=argparse.ArgumentParser(description="Prepare submission scripts for sge, pbs, or slurm cluster management systems.")

swd=os.path.dirname(os.path.realpath(__file__))

ap.add_argument("--regions", help="Full path to regions.", required=True)
ap.add_argument("--params", help="Full path to parameter file.", required=True)
ap.add_argument("--runTrio", help="Prepare commands for trio assembly.", default=False, action='store_true')
ap.add_argument("--grid", help="Grid type", choices=["sge", "pbs", "slurm"], required=True)
ap.add_argument("--base", help="base name for grid commands", default="psv_grid")
ap.add_argument("--conc", help="Number of concurrent jobs", default="50")
ap.add_argument("--grouped", help="(SLURM only) Some cluster configurations do not allow very many array jobs. Normally execution is one region assembled per job, which enables fine-scale job management at the cost of high load on the CMS. This instead makes one job run many (~1000) regions per job. Since each job takes up to 5 minutes, this can take up to 3-4 days per job.", default=0, required=False, type=int)
ap.add_argument("--config", help="Extra configuration parameters for job submission, for example specification of a specific queue, etc", default="")

args=ap.parse_args()

regionsFile = open(args.regions)
numRegions = len(regionsFile.readlines())
assemblyScript=swd+"/../RunTiledAssembly.sh"
if args.runTrio:
    assemblyScript=swd+"/../RunTrioTiledAssemblyOnRegions.sh"


if args.grid == "slurm":
    gridFileName=args.base+".run.sh"
    gf=open(gridFileName,'w')
    gf.write("#!/usr/bin/env bash\n")
    nBins = 0
    if args.grouped > 0:
        regionsPerJob = numRegions // args.grouped +  1
        gf.write("for i in {{1..{}}}; do\n".format(regionsPerJob))
        gf.write("   index=$(($i*{} + $SLURM_ARRAY_TASK_ID))\n".format(args.grouped))
        gf.write("   " + assemblyScript +" `awk \"NR == $index\" " + args.regions + " ` " + args.params + "\n")       
        gf.write("done\n")
        gf.close()

        submit=args.base+".submit.sh"
        sf=open(submit,'w')
        sf.write("#!/usr/bin/env bash\n")
        sf.write("sbatch --cpus-per-task=4 --mem=4G --time=72:00:00 --array=0-{}%{} {} {}\n".format(args.grouped-1, args.conc, gridFileName, args.config))
        sf.close()

    else:
        gf.write(assemblyScript +" `awk \"NR == $SLURM_ARRAY_TASK_ID\" " + args.regions + " ` " + args.params + "\n")
        gf.close()

        submit=args.base+".submit.sh"
        sf=open(submit,'w')
        sf.write("#!/usr/bin/env bash\n")
        sf.write("sbatch --cpus-per-task=4 --mem=4G --time=1:00:00 --array=1-{}%{} {} {}\n".format(numRegions, args.conc, gridFileName, args.config))
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
    
    
def GetMaxTasksPBS():
    sys.stderr.write("""CAUTION! Setting max tasks to 10,000 although this may not 
be compatible with your PBS system.""")
    return 10000

if args.grid == "pbs":
    maxJobsPerFile=GetMaxTasksPBS()

    nFiles = int(math.ceil(float(numRegions)/maxJobsPerFile))
    startJob=1
    allJobs=[]
    for i in range(0,nFiles):
        gridFileName=args.base + "." + str(i) + ".run.sh"
        gf=open(gridFileName,'w')
        gf.write("#!/usr/bin/env bash\n")
        endJob=min(numRegions, startJob+maxJobsPerFile-1)
        gf.write("#PBS -t {}-{}%{}\n".format(startJob, endJob, args.conc))
        gf.write("#PBS -l walltime=30:00:00,nodes=1:ppn=4,mem=4gb\n")
        gf.write(assemblyScript +" `awk \"NR == $PBS_ARRAYID\" " + args.regions + " ` " + args.params + "\n")
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
