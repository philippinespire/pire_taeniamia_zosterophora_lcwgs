#!/bin/bash -l

#SBATCH --job-name=angsd_sfs
#SBATCH -o angsd_sfs-%j.out 
#SBATCH -p main 
#SBATCH -n 1  
#SBATCH -N 1  
#SBATCH --cpus-per-task=8 

enable_lmod
module load container_env angsd
export SINGULARITY_BIND=/home/e1garcia

OUTDIR=$1


cd ${OUTDIR}

crun.angsd realSFS abol_sites_notrans.saf.idx -P 8 -fold 1 > abol_notrans.sfs
