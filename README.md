# pire_taeniamia_zosterophora_lcwgs
Repository for bioinformatics and analysis of historical and contemporary low-coverage whole genome sequencing data of Taeniamia zosterophora collected in the Philippines.

# Pre-processing
Following the pire_fq_gz_processing repo

Files were already downloaded, and the decode file looks fine, so moving to Step 5: renaming

# 5. Perform a renaming dry-run
   Run renameFQGZ_keeplane.bash to view the original and new file names and create tsv files to store the original and new file naming conventions.
   Need to use the keeplane script because each sample was sequenced across multiple lanes.
 
log into a compute node interactively so this goes faster
salloc

once you have the compute node, procede
bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/renameFQGZ_keeplane.bash Tzo-LCWGS-FullSeq-SequenceNameDecode.tsv

Dry run looks okay, so will proceed to renaming for real.

# 6. Rename the files for real
bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/renameFQGZ_keeplane.bash Tzo-LCWGS-FullSeq-SequenceNameDecode.tsv rename

#you will need to say y 2X

# 7. Check quality of the data with fastqc

Run fastqc and multiqc with the Multi_FASTQC.sh script 

sbatch /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/Multi_FASTQC.sh "/archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/fq_raw" "fqc_raw_report" "fq.gz"

Initial quality screen: looks good overall, GC content seems reasonable, percent duplication seems workable
Potential issues:  
  * % duplication - 
	* Alb: 5-20%, Contemp: 5-20%
  * GC content - 
	* Alb: 40-55%, Contemp: 40-55%
  * number of reads - 
	* Alb: 0.5-20 mil, Contemp: 0.5-20 mil

# 8. First trim

Execute runFASTP_1st_trim.sbatch (0.5-3 hours run time)

sbatch /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/runFASTP_1st_trim.sbatch fq_raw fq_fp1

After 1st trim: overall looks okay, only concern is a high % adapter in some samples
Potential issues:  
  * % duplication - 
    * Alb: 1-14%, Contemp: 1-14%
  * GC content -
    * Alb: 36-47%, Contemp: 40-45 %
  * passing filter - 
    * Alb: 91-98%, Contemp: 82-98%
  * % adapter - 
    * Alb: 40-98%, Contemp: 3-82%
  * number of reads - 
    * Alb: 11-4.5K mil, Contemp: 100-8K mil

# 9. Remove duplicates with clumpify

Execute runCLUMPIFY_r1r2_array.bash (0.5-3 hours run time)**

runCLUMPIFY_r1r2_array.bash is a bash script that executes several sbatch jobs to de-duplicate and clumpify your fq.gz files. It does two things:

    Removes duplicate reads.
    Re-orders each fq.gz file so that similar sequences (reads) appear closer together. This helps with file compression and speeds up downstream steps.

You will need to specify the number of nodes you wish to allocate your jobs to. The max # of nodes to use at once should not exceed the number of pairs of r1-r2 files to be processed. (Ex: If you have 3 pairs of r1-r2 files, you should only use 3 nodes at most.) If you have many sets of files (likely to occur if you are processing capture data), you might also limit the nodes to the current number of idle nodes to avoid waiting on the queue (run sinfo to find out # of nodes idle in the main partition)

bash /archive/carpenterlab/pire/pire_fq_gz_processing/runCLUMPIFY_r1r2_array.bash fq_fp1 fq_fp1_clmp /scratch/kfitz012 40

Generate metadata on deduplicated FASTQ files

Once CLUMPIFY has finished running and there are no issues, run runMULTIQC.sbatch to get the MultiQC output.

sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/Multi_FASTQC.sh "fq_fp1_clmp" "fqc_clmp_report"  "fq.gz"

Results look good, duplication rate is low now

# 10. Second trim

Run runFASTP_2_cssl.sbatch script

/home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runFASTP_2_cssl.sbatch fq_fp1_clmp fq_fp1_clmp_fp2 33

Modified script to a front trim of 30bp so that fastp would run.

MultiQC output looks good, seems to have removed most duplication. GC content has decreased a little, but is still around 40%, which seems reasonable. 
Potential issues:
% duplication -
Alb: 0-2%%, Contemp: 0-2%

GC content -
Alb: 38-45%, Contemp: 38-45%

passing filter -
Alb: 33-95%, Contemp: 75-98%

% adapter -
Alb: 3-45%, Contemp: 3-40%

number of reads -
Alb: 0.1-25 mil, Contemp: 0.1-50 mil

Overall seems good to move onto decontamination.

# 11. Decontamination

Run runFQSCRN_6.bash script

bash

fqScrnPATH=/home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runFQSCRN_6.bash

indir=fq_fp1_clmp_fp2

outdir=/scratch/kfitz012/fq_fp1_clmp_fp2_fqscrn

nodes=12 #Wahab is saying that 12 nodes @ 40 threads per node is the max CPU per user

bash $fqScrnPATH $indir $outdir $nodes


