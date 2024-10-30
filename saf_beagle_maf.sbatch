#!/bin/bash -l

#SBATCH --job-name=angsd_saf
#SBATCH -o angsd_saf-%j.out 
#SBATCH -p main 
#SBATCH -n 1  
#SBATCH -N 1  
#SBATCH --cpus-per-task=8 

enable_lmod
module load container_env angsd
export SINGULARITY_BIND=/home/e1garcia

OUTDIR=$1


cd ${OUTDIR}

crun angsd \
        -b bam_list_abol.txt \
        -anc /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/reference/reference.denovoSSL.Tzo20k.fasta \
        -ref /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/reference/reference.denovoSSL.Tzo20k.fasta \
        -out abol_sites_notrans \
      -doSaf 1 -noTrans 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1  \
        -P 8 \
	-sites global_snp_list_neutral.txt -rf global_snp_list_depth1_15.chrs
