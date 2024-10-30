#!/bin/bash -l

#SBATCH --job-name=angsd_fst
#SBATCH -o angsd_fst-%j.out 
#SBATCH -p main 
#SBATCH -n 1  
#SBATCH -N 1  
#SBATCH --cpus-per-task=8 

enable_lmod
module load container_env angsd
export SINGULARITY_BIND=/home/e1garcia


## This script is used to get pairwise Fst estimates from angsd for each population / group pair

#This script generates a number of files, so creating a new directory and copying the .saf.idx files into it.

SAFDIR=${1:-/archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/angsd/bam/fst 
THREADS=${2:-8} # Number of parallel threads to use, default value is 8.
EXTRA_ARG=${3:-''} # Extra arguments for the SFS estimation step, default value is ''
POP1=${4:-abol}
POP2=${5:-amta}

cd $SAFDIR

 # Generate the 2dSFS to be used as a prior for Fst estimation (and individual plots)
  crun.angsd realSFS $POP1'.saf.idx' $POP2'.saf.idx' -P $THREADS $EXTRA_ARG > $POP1'_'$POP2'.2dSFS'

 # Estimating Fst in angsd
 crun.angsd realSFS fst index  $POP1'.saf.idx' $POP2'.saf.idx' -sfs $POP1'_'$POP2'.2dSFS' -fstout $POP1'_'$POP2'.alpha_beta'
 crun.angsd realSFS fst print $POP1'_'$POP2'.alpha_beta.fst.idx' > $POP1'_'$POP2'.alpha_beta.txt'
  awk '{ print $0 "\t" $3 / $4 }' $POP1'_'$POP2'.alpha_beta.txt' > $POP1'_'$POP2'.fst'

 # Estimating average Fst in angsd
 crun.angsd realSFS fst stats $POP1'_'$POP2'.alpha_beta.fst.idx' > $POP1'_'$POP2'.average_fst.txt'
