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
OutFilePrefix=$(basename "$OutPrefix")

ResultsDir="${OutDir}/QTL.Results"

mkdir -p $OutDir
mkdir -p $ResultsDir

if [ $chr == all ] 
then 
  chr=($(seq 1 1 22))   # seq FIRST STEP LAST
else
  IFS=',' read -r -a chr <<< "$chr"
fi

for i in ${chr[@]} 
do
  echo "***************************************************************************************"
  echo "                               QTL analysis chromosome $i ..."
  echo "***************************************************************************************"
  Rscript ${ScriptDir}/Main.R "${ResultsDir}/${OutFilePrefix}" $i $Dist $CisPvalue $TransPvalue $TransCrossChr "$Overwrite"
done



