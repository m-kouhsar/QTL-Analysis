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

OutDir=$(dirname "$OutPrefix")
mkdir -p $OutDir

if [ $chr == all ] 
then 
  chr=($(seq 1 1 22))   # seq FIRST STEP LAST
else
  IFS=',' read -r -a chr <<< "$chr"
fi

for i in ${chr[@]} 
do
  echo "---------------------------------------------------------\n"
  echo "Working on chromosome $i ..."
  Rscript ${ScriptDir}/Main.R $OutPrefix $i $Dist $CisPvalue $TransPvalue $TransCrossChr 
done



