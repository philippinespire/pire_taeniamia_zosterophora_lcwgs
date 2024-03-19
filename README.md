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

Run runFQSCRN_6.bash script (had to adjust script to run with 20 threads/node because otherwise only 12 jobs could run at a time)

bash

fqScrnPATH=/archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/runFQSCRN_6.bash

indir=fq_fp1_clmp_fp2

outdir=/scratch/kfitz012/fq_fp1_clmp_fp2_fqscrn

nodes=24 #Wahab is saying that 24 nodes @ 20 threads per node is the max CPU per user

bash $fqScrnPATH $indir $outdir $nodes


MultiQC:

sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runMULTIQC.sbatch fq_fp1_clmp_fp2_fqscrn fastq_screen_report

MultiQC output looks good. Contamination was low (<5% of reads) for both historical and contemporary samples.
 
Potential issues:
one hit, one genome, no ID -
Alb: 95-97%, Contemp: 95-97%

no one hit, one genome to any potential contaminators (bacteria, virus, human, etc) -
Alb: 3-5%, Contemp: 3-5%

Re-pair files 

sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runREPAIR.sbatch /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/fq_fp1_clmp_fp2_fqscrn/fq_fp1_clmp_fp2_fqscrn /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/fq_fp1_clmp_fp2_fqscrn_rprd 5

Check that the files re-paired properly

bash

SCRIPT=/home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/validateFQPE.sbatch

DIR=fq_fp1_clmp_fp2_fqscrn_rprd

fqPATTERN="*fq.gz"

sbatch $SCRIPT $DIR $fqPATTERN

Output looks good. The only issue is 7 files did not have any reads, but this lines up with the earlier MultiQC outputs where some sequences did not have any reads other than the adapter.
Since samples were sequenced across lanes, no individuals have been lost at this point.

MultiQC with the re-paired files:

#sbatch Multi_FASTQC.sh "<indir>" "<output report name>" "<file extension>"
sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/Multi_FASTQC.sh "/archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/fq_fp1_clmp_fp2_fqscrn_rprd" "fqc_rprd_report" "fq.gz"

MultiQC results look good. Duplication remains low and GC content seems appropriate for a fish species. Some reads are shorter than others but that is to be expected with historical DNA.
 Sequence files are ready for mapping.

Potential issues:
% duplication -
Alb: 0-4%, Contemp: 0-4%

GC content -
Alb: 36-43%, Contemp: 36-43%

number of reads -
Alb: 0.5-14 mil, Contemp: 0.5-15 mil


# Moving to the pire_lcwgs_data_processing pipeline for mapping and down-stream analysis

# 2 & 3 Get and curate the reference genome

Using the same reference genome as used for Tzo CSSL work. From the Tzo CSSL README: Found the best genome by running wrangleData.R, sorted tibble by busco single copy complete,
quast n50, and filtered by species in Rstudio. The best genome to map for Tzo is Tzo_scaffolds_TzC0402G_contam_R1R2_noIsolate.fasta in 
/home/e1garcia/shotgun_PIRE/pire_ssl_data_processing/taeniamia_zosterophora/probe_design.

Moving the reference genome to my working directory:

mkdir refGenome
cd refGenome
cp /home/e1garcia/shotgun_PIRE/pire_ssl_data_processing/taeniamia_zosterophora/probe_design/Tzo_scaffolds_TzC0402G_contam_R1R2_noIsolate.fasta /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/refGenome

#the destination reference fasta should be named as follows: reference.<assembly type>.<unique assembly info>.fasta
#<assembly type> is `ssl` for denovo assembled shotgun library or `rad` for denovo assembled rad library
#this naming is a little messy, but it makes the ref 100% tracable back to the source
#it is critical not to use `_` in name of reference for compatibility with ddocent and freebayes

mv Tzo_scaffolds_TzC0402G_contam_R1R2_noIsolate.fasta ./reference.ssl.Tzo-C-0402G-R1R2-contam-noisolate.fasta


# 4 Map reads to reference genome 

#sbatch mkBAMwgs.sbatch InFilePattern RefGenomeFile OutDir

sbatch /archive/carpenterlab/pire/pire_lcwgs_data_processing/scripts/mkBAM.sbatch "/archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/fq_fp1_clmp_fp2_fqscrn_rprd/*.fq.gz" /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/refGenome/reference.ssl.Tzo-C-0402G-R1R2-contam-noisolate.fasta /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/mkBAM


