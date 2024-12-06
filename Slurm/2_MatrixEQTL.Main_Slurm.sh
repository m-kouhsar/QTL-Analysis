#!/bin/bash
#SBATCH -A Research_Project1 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=00:05:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

source $1

if [ $chr == all ]
then
  array=($(seq 1 1 22))  ## seq FIRST STEP LAST
  printf -v chr '%s,' "${array[@]}"
  chr="${chr%,}"
fi


sbatch --array=$chr ${ScriptDir}/Main.sh ${ScriptDir} $OutPrefix $Dist $CisPvalue $TransPvalue $TransCrossChr $Overwrite



