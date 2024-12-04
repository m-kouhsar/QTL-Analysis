#!/bin/bash

source $1

Overwrite="${Overwrite#[[:space:]]}"
Overwrite="${Overwrite%[[:space:]]}"
Overwrite="$(echo "$Overwrite" | tr '[:upper:]' '[:lower:]')"

chr="${chr#[[:space:]]}"
chr="${chr%[[:space:]]}"
chr="$(echo "$chr" | tr '[:upper:]' '[:lower:]')"

OutDir=$(dirname "$OutPrefix")
mkdir -p $OutDir

IFS=',' read -r -a array <<< "$chr"

if [ $chr == all ]
then
  array=($(seq 1 1 22))  ## seq FIRST STEP LAST
fi

for i in "${array[@]}"
do
  if [ ! -f ${OutPrefix}.chr${i}.raw ] || [ $Overwrite = "yes" ]
  then
    echo "Converting binary format chromosome $i..."
    plink --bfile ${GenotypeBinaryPrefix} --recodeA --chr $i --out ${OutPrefix}.chr${i}
    echo "#########################################################################################"
  fi
done

if [ ! -f ${OutPrefix}.eigenvec ] || [ $Overwrite = "yes" ]
then
  echo "Calculating genotype data PCs..."
  gcta64 --bfile ${GenotypeBinaryPrefix} --make-grm-bin --out ${OutPrefix} --thread-num 16
  gcta64 --grm ${OutPrefix} --pca --out ${OutPrefix}
fi


Rscript ${ScriptDir}/PrepareData.R ${GenotypeBinaryPrefix}.fam ${OutPrefix}.eigenvec $ExpressionFile $PhenotypeFile $GeneLocationFile "$FactCovar" "$NumCovar" "$OutPrefix" $chr $Overwrite

