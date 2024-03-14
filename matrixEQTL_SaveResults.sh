#!/bin/sh
#!/bin/bash
#SBATCH -A Research_Project1 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

### print start date and time
echo Job started on:
date -u

set -e
####### 

### NOTE: Do not store confidenial information in this file use the config file

######
InDir=./Results/UKBBN
FilePrefix=UKBBN
Dist=250000
cis_pval=0.05
trans_pval=0.05
trans_cross_chr=no
save_csv_cis=yes
save_csv_trans=yes

SCRIPTDIR=./R

module load R

Rscript ${SCRIPTDIR}/matrixQTL_eQTL_SaveResults.R $InDir $FilePrefix $cis_pval $trans_pval $Dist $trans_cross_chr $save_csv_cis $save_csv_trans

echo Job finished:
date -u
