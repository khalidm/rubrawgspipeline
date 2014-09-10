#!/bin/bash
# Created by the VLSCI job script generator for SLURM on x86
# Fri Jan 24 2014 15:05:28 GMT+1100 (EST)

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
       echo "You need to submit your job to the queuing system with sbatch"
          exit 1
      fi

# Partition for the job:
#SBATCH -p main

# Maximum number of CPU cores used by the job:
#SBATCH --ntasks=1

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-4:0:00

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL

# Use this email address:
#SBATCH --mail-user=khalid.mahmood@unimelb.edu.au

# Run the job from the directory where it was launched (default):
# The job command(s):
cat $in | gzip -9 > $out

