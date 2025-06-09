#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=10:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=eQTM.%j.out

############################################################################

# DNAm_file: A tsv file contains normalized methylation values. Rows must be CpG IDs and columns must be Sample IDs
# DNA_annot_file: CpG annotation file. It must contains the following columns:
#                         id: CpG IDs
#                         chr: Chromosome
#                         start: Start position
#                         end: End position
# RNAe_file: A tsv file contains normalized expression values. Rows must be RNA/Gene IDs and columns must be Sample IDs
# RNA_annot_file: RNA/Gene annotation file. It must contains the following columns:
#                         id: RNA/Gene IDs
#                         chr: Chromosome
#                         start: Start position
#                         end: End position
# Pheno_file: A tsv file contains phenotype information about samples
# covars_fact: Categorical covariates (eg. Sex) for adding to the lm model (must be represented in column names in Pheno_file). 
#              The first categorical variable will be used for coloring the Plots 
# covars_num: Numerical covariates (eg. Age) for adding to the lm model (must be represented in column names in Pheno_file) 
# d_cis: Distance from start and end of each RNA/Gene to search for CpGs
# PCs: Numebr of Principal Components (in both methylation and expression data) need to be added to the lm model
# OutPrefix: Result file prefix

############################################################################

DNAm_file="/lustre/projects/Research_Project-191391/Morteza/kallisto/eQTM/BDR.AD.C.DNAm.tsv"
DNA_annot_file="/lustre/projects/Research_Project-191391/Morteza/kallisto/eQTM/DNAm.EPIC.Annotation.hg38.tsv"
RNAe_file="/lustre/projects/Research_Project-191391/Morteza/kallisto/eQTM/BDR.AD.C.RNAe.GRCh38.111.tpm.tsv"
RNA_annot_file="/lustre/projects/Research_Project-191391/Morteza/kallisto/eQTM/RNA.Gene.Annot.GRCh38.111.tsv"
Pheno_file="/lustre/projects/Research_Project-191391/Morteza/kallisto/eQTM/BDR.AD.C.Phenotype.csv"
covars_fact="Phenotype,Gender"
covars_num="Age,RIN"
d_cis=25000
PCs=10
OutPrefix="/lustre/projects/Research_Project-191391/Morteza/kallisto/eQTM/Results/BDR"

ScriptDir="/lustre/projects/Research_Project-191391/Morteza/kallisto/eQTM/Scripts"


Rscript ${ScriptDir}/eQTM.R $DNAm_file $DNA_annot_file  $RNAe_file  $RNA_annot_file $Pheno_file $covars_fact $covars_num $d_cis $PCs $OutPrefix 

