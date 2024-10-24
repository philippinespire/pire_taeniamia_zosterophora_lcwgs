#!/bin/bash -l

#SBATCH --job-name=ATLAS_recal_theta
#SBATCH -o ATLAS_recal_theta-%j.out 
#SBATCH -p main 
#SBATCH -n 1  
#SBATCH -N 1  
#SBATCH --cpus-per-task=20 

enable_lmod
module load container_env bwa
module load container_env samtools
module load container_env atlas/20240503
export SINGULARITY_BIND=/home/e1garcia

BAMDIR=$1
REFLOC=$2
OUTDIR=$3

BAMPATTERN=*merged.rmdup.merged.fltrd.bam

all_samples=$(ls $BAMDIR/$BAMPATTERN | sed -e 's/\.merged\.rmdup\.merged\.fltrd\.bam//' -e 's/.*\///g')
all_samples=($all_samples)

sample_name=${all_samples[${SLURM_ARRAY_TASK_ID}]}
echo ${sample_name}

cd ${OUTDIR}

cp ${BAMDIR}/${sample_name}.merged.rmdup.merged.fltrd.bam .

crun.samtools samtools index ${sample_name}.merged.rmdup.merged.fltrd.bam

crun.atlas atlas --task estimateErrors --bam ${sample_name}.merged.rmdup.merged.fltrd.bam --fasta ${REFLOC} --regions gerp15.merge.nosingle.bed,gerp1.merge.nosingle.bed,gerp5.merge.nosingle.bed --readUpToDepth 10
