#!/bin/bash -l

#SBATCH --job-name=ATLAS_convertvcf_recal
#SBATCH -o ATLAS_convertvcf_recal-%j.out 
#SBATCH -p main 
#SBATCH -n 1  
#SBATCH -N 1  
#SBATCH --cpus-per-task=40 

enable_lmod
module load container_env atlas/20240503
export SINGULARITY_BIND=/home/e1garcia

OUTDIR=$1


cd ${OUTDIR}


crun.atlas atlas --task convertVCF --vcf majorMinor_0.1_5.vcf.gz --format beagle --minMAF 0.001 --out atlas_0.1_5_beagle --logFile beagle_0.1_5.log
