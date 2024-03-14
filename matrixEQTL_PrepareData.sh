#!/bin/bash
#SBATCH -A Research_Project1 # research project to submit under.
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
echo ""
####### 

SCRIPTDIR=./R
InDir=./Raw
OutDir=./Results
FilePrefix=UKBBN
covar_fact=""
covar_num=""
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

mkdir -p $OutDir/${FilePrefix}

IFS=',' read -r -a array <<< "$chr"

if [ $chr == all ]
then
  array=($(seq 1 1 22))  ## seq FIRST STEP LAST
fi

for i in "${array[@]}"
do
  if [ ! -f ${OutDir}/${FilePrefix}/${FilePrefix}.chr${i}.raw ]
  then
    echo "Running chr ${i}..."
    plink --bfile ${InDir}/${FilePrefix} --recodeA --chr $i --out ${OutDir}/${FilePrefix}/${FilePrefix}.chr${i}
    echo "#########################################################################################"
  fi
done

if [ ! -f ${OutDir}/${FilePrefix}/${FilePrefix}.eigenvec ]
then
  gcta64 --bfile ${InDir}/${FilePrefix} --make-grm-bin --out ${OutDir}/${FilePrefix}/${FilePrefix} --thread-num 16
  gcta64 --grm ${OutDir}/${FilePrefix}/${FilePrefix} --pca --out ${OutDir}/${FilePrefix}/${FilePrefix}
fi

Rscript ${SCRIPTDIR}/matrixQTL_eQTL_PrepareData.r $InDir $OutDir/${FilePrefix} "$covar_fact" "$covar_num" "$FilePrefix" $chr


echo Job finished:
date -u
