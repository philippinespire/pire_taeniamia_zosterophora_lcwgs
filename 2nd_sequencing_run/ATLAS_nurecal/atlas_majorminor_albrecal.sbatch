#!/bin/bash -l

#SBATCH --job-name=ATLAS_majorminor_recal
#SBATCH -o ATLAS_majorminor_recal-%j.out 
#SBATCH -p main 
#SBATCH -n 1  
#SBATCH -N 1  
#SBATCH --cpus-per-task=40 

enable_lmod
module load container_env atlas/20240503
export SINGULARITY_BIND=/home/e1garcia

OUTDIR=$1


cd ${OUTDIR}


crun.atlas atlas --task majorMinor --glf glf_0.1_5.txt --minMAF 0.001 --out majorMinor_0.1_5 --logFile majorMinor_0.1_5.log
