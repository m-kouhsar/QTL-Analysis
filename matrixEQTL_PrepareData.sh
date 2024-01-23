#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=2:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

### print start date and time 
echo Job started on:
date -u
###

SCRIPTDIR=./R
InDir=./Inputs
OutDir=./Outputs
FilePrefix=UKBBN
covar_fact=Gender,Brain.Bank,RINcat,Plate
covar_num=Age
chr=all  ##for multiple chr use chr=2,3,5,8,19 or chr=all

echo "Scritp directory: "$SCRIPTDIR
echo "Input directory: "$InDir
echo "Output directory: "$OutDir
echo "Input and Output files prefix: "$FilePrefix
echo "Factor covariates: "$covar_fact
echo "Numeric covariates: "$covar_num
echo "Chromosome: "$chr
echo "############################################################"
echo ""

mkdir -p $OutDir

module load R

IFS=',' read -r -a array <<< "$chr"

if [ $chr == all ]
then
  array=($(seq 1 1 22))  ## seq FIRST STEP LAST
fi

for i in "${array[@]}"
do
  if [ ! -f ${OutDir}/${FilePrefix}.chr${i}.raw ]
  then
    echo "Running chr ${i}..."
    plink --bfile ${InDir}/${FilePrefix} --recodeA --chr $i --out ${OutDir}/${FilePrefix}.chr${i}
    echo "#########################################################################################"
  fi
done

if [ ! -f ${OutDir}/${FilePrefix}.eigenvec ]
then
  gcta64 --bfile ${InDir}/${FilePrefix} --make-grm-bin --out ${OutDir}/${FilePrefix} --thread-num 16
  gcta64 --grm ${OutDir}/${FilePrefix} --pca --out ${OutDir}/${FilePrefix}
fi

Rscript ${SCRIPTDIR}/matrixQTL_eQTL_PrepareData.r $InDir $OutDir $covar_fact $covar_num $FilePrefix $chr


echo Job finished:
date -u
