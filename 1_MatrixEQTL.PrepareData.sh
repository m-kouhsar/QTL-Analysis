#!/bin/bash

source $1

Overwrite="${Overwrite#[[:space:]]}"
Overwrite="${Overwrite%[[:space:]]}"
Overwrite="$(echo "$Overwrite" | tr '[:upper:]' '[:lower:]')"

chr="${chr#[[:space:]]}"
chr="${chr%[[:space:]]}"
chr="$(echo "$chr" | tr '[:upper:]' '[:lower:]')"

OutDir=$(dirname "$OutPrefix")
OutFilePrefix=$(basename "$OutPrefix")
PlinkDir="${OutDir}/QTL.PlinkData"

mkdir -p $PlinkDir

IFS=',' read -r -a array <<< "$chr"

if [ $chr == all ]
then
  array=($(seq 1 1 22))  ## seq FIRST STEP LAST
fi

for i in "${array[@]}"
do
  if [ ! -f ${PlinkDir}/${OutFilePrefix}.chr${i}.raw ] || [ $Overwrite = "yes" ]
  then
    echo "Converting binary format chromosome $i..."
    plink --bfile ${GenotypeBinaryPrefix} --recodeA --chr $i --out ${PlinkDir}/${OutFilePrefix}.chr${i}
    echo "#########################################################################################"
  else
    echo "Formatted genotype data for chromosome $i already exist."
    echo "#########################################################################################"
  fi
done

if [ ! -f ${PlinkDir}/${OutFilePrefix}.eigenvec ] || [ $Overwrite = "yes" ]
then
  echo "Calculating genotype data PCs..."
  gcta64 --bfile ${GenotypeBinaryPrefix} --make-grm-bin --out ${PlinkDir}/${OutFilePrefix} --thread-num 16
  gcta64 --grm ${PlinkDir}/${OutFilePrefix} --pca --out ${PlinkDir}/${OutFilePrefix}
else
    echo "Genotype data PCs already exists."
    echo "#########################################################################################"
fi


Rscript ${ScriptDir}/PrepareData.R ${GenotypeBinaryPrefix}.fam ${PlinkDir}/${OutFilePrefix}.eigenvec $ExpressionFile $PhenotypeFile $GeneLocationFile "$FactCovar" "$NumCovar" "$OutPrefix" "$chr" "$Overwrite"

