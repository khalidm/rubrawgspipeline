#!/bin/bash
# Created by the VLSCI job script generator for SLURM on x86
# Tue Jan 21 2014 15:57:44 GMT+1100 (EST)

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Partition for the job:
#SBATCH -p main


# JOBNAME=$NAME
# The name of the job:
# #SBATCH --job-name=B2F

# Maximum number of CPU cores used by the job:
# #SBATCH --ntasks=1

# The amount of memory in megabytes per process in the job:
# #SBATCH --mem-per-cpu=65536
#SBATCH --mem=65536

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-36:0:00

# SMP job: want exclusive use of 1 node
#SBATCH --nodes=1
#SBATCH --exclusive

# #SBATCH --output=logfiles/%j.out
# #SBATCH --error=logfiles/%j.err

# Run the job from the directory where it was launched (default):
# The modules to load:
module load java/1.7.0_25
module load picard/1.69

# The job command(s):
# mem=$1
# shift

java -d64 -Xmx24g -Djava.io.tmpdir=$TMPDIR -jar $PICARD_HOME/lib/SamToFastq.jar TMP_DIR=$TMPDIR INPUT=$in OUTPUT_DIR=$out OUTPUT_PER_RG=TRUE

# /vlsci/VR0182/shared/km/pipeline/PicardSamToFastq $INPUT

