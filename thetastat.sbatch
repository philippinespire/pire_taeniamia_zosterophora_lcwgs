#!/bin/bash -l

#SBATCH --job-name=angsd_theta
#SBATCH -o angsd_theta-%j.out 
#SBATCH -p main 
#SBATCH -n 1  
#SBATCH -N 1  
#SBATCH --cpus-per-task=8

enable_lmod
module load container_env angsd
export SINGULARITY_BIND=/home/e1garcia

OUTDIR=$1


cd ${OUTDIR}


crun.angsd thetaStat do_stat abol_notrans.thetas.idx
