#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

source $1

OutDir=$(dirname "$OutPrefix")
mkdir -p $OutDir

IFS=',' read -r -a array <<< "$chr"

if [ $chr == all ]
then
  array=($(seq 1 1 22))  ## seq FIRST STEP LAST
fi

for i in "${array[@]}"
do
  if [ ! -f ${OutPrefix}.chr${i}.raw ]
  then
    echo "Running chr ${i}..."
    plink --bfile ${GenotypeBinaryPrefix} --recodeA --chr $i --out ${OutPrefix}.chr${i}
    echo "#########################################################################################"
  fi
done

if [ ! -f ${OutPrefix}.eigenvec ]
then
  gcta64 --bfile ${GenotypeBinaryPrefix} --make-grm-bin --out ${OutPrefix} --thread-num 16
  gcta64 --grm ${OutPrefix} --pca --out ${OutPrefix}
fi


Rscript ${ScriptDir}/PrepareData.R ${GenotypeBinaryPrefix}.fam ${OutPrefix}.eigenvec $ExpressionFile $PhenotypeFile $GeneLocationFile "$FactCovar" "$NumCovar" "$OutPrefix" $chr

