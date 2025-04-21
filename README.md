# pire_taeniamia_zosterophora_lcwgs
Repository for bioinformatics and analysis of historical and contemporary low-coverage whole genome sequencing data of Taeniamia zosterophora collected in the Philippines.

This repository provides the code for the paper: Preservation of genetic diversity and selection over a century in a coral reef fish (Taeniamia zosteorphora) in the Philippines.

As genomic data is too large to be stored on GitHub, sequence data will be uploaded to the NCBI SRA and metadata to GEOME-DB.

Metadata Links:
Malampaya: https://n2t.net/ark:/21547/FtD2 (Historical) and https://n2t.net/ark:/21547/FtS2 (Contemporary)  
Mantatao Island: https://n2t.net/ark:/21547/FtY2 (Historical) and https://n2t.net/ark:/21547/Fte2 (Contemporary)

Sequence Data Links and NCBI BioProject Numbers:

Malampaya Historical: PRJNA1185294 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1185294)  
Malampaya Contemporary: PRJNA1185280 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1185280)  
Mantatao Island Historical: PRJNA1185271 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1185271)  
Mantatao Island Contemporary: PRJNA1185264 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1185264)

The genome assembly used for mapping is available on FigShare: 10.6084/m9.figshare.27682767.

The README is separated into 4 main sections:
 1) Pre-processing of sequence data following the pire_fq_gz_processing repository (https://github.com/philippinespire/pire_fq_gz_processing)
 2) Mapping using the GenErode pipeline (generating rescaled and realigned bam files for use in ANGSD)
 3) Downstream analysis using ANGSD
 4) Statistical analysis in R

Intermediate outputs for pre-processing data steps are found in the 1st_sequencing_run and 2nd_sequencing_run directories (https://github.com/philippinespire/pire_taeniamia_zosterophora_lcwgs/tree/main/1st_sequencing_run and https://github.com/philippinespire/pire_taeniamia_zosterophora_lcwgs/tree/main/2nd_sequencing_run).

Config files and scripts for mapping using the GenErode pipeline are found in the 2nd_sequencing_run/GenErode_20k_2runs_2 directory (https://github.com/philippinespire/pire_taeniamia_zosterophora_lcwgs/tree/main/2nd_sequencing_run/GenErode_20k_2runs_2).

Downstream analysis using ANGSD and Statistical analysis in R used scripts that can be found in the Scripts directory (https://github.com/philippinespire/pire_taeniamia_zosterophora_lcwgs/tree/main/Scripts), and select intermediate files can be found in the Files directory (https://github.com/philippinespire/pire_taeniamia_zosterophora_lcwgs/tree/main/Files). 

Pre-processing, mapping, and downstream analysis were run using Wahab, Old Dominion University's high performance computing system.  
Statistical analysis was run using R v.4.4.1 on a MacBook Pro (Apple M3 Pro Chip with 36 GB memory)

# 1. Pre-processing
Following the pire_fq_gz_processing repo (https://github.com/philippinespire/pire_fq_gz_processing)

Files were already downloaded, and the decode file looks fine, so moving to Step 5: renaming

## Perform a renaming dry-run
Run renameFQGZ_keeplane.bash to view the original and new file names and create tsv files to store the original and new file naming conventions.
Need to use the keeplane script because each sample was sequenced across multiple lanes.
 
log into a compute node interactively so this goes faster
salloc

once you have the compute node, procede
bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/renameFQGZ_keeplane.bash Tzo-LCWGS-FullSeq-SequenceNameDecode.tsv

Dry run looks okay, so will proceed to renaming for real.

## Rename the files for real
bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/renameFQGZ_keeplane.bash Tzo-LCWGS-FullSeq-SequenceNameDecode.tsv rename


## Check quality of the data with fastqc

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

## First trim

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

## Remove duplicates with clumpify

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

## Second trim

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

## Decontamination

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

sbatch Multi_FASTQC.sh "<indir>" "<output report name>" "<file extension>"
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


## Moving to the pire_lcwgs_data_processing pipeline (https://github.com/philippinespire/pire_lcwgs_data_processing) to curate the reference genome assembly for mapping with the GenErode pipeline

## Get and curate the reference genome

Using the same reference genome as used for Tzo CSSL work (https://github.com/philippinespire/pire_taeniamia_zosterophora_cssl). From the Tzo CSSL README: Found the best genome by running wrangleData.R, sorted tibble by busco single copy complete,
quast n50, and filtered by species in Rstudio. The best genome to map for Tzo is Tzo_scaffolds_TzC0402G_contam_R1R2_noIsolate.fasta in 
/home/e1garcia/shotgun_PIRE/pire_ssl_data_processing/taeniamia_zosterophora/probe_design.

If not using Wahab, the genome assembly can be found on FigShare: 10.6084/m9.figshare.27682767.

Moving the reference genome to my working directory:

mkdir refGenome
cd refGenome
cp /home/e1garcia/shotgun_PIRE/pire_ssl_data_processing/taeniamia_zosterophora/probe_design/Tzo_scaffolds_TzC0402G_contam_R1R2_noIsolate.fasta /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/refGenome

#the destination reference fasta should be named as follows: reference.<assembly type>.<unique assembly info>.fasta
#<assembly type> is `ssl` for denovo assembled shotgun library or `rad` for denovo assembled rad library
#this naming is a little messy, but it makes the ref 100% tracable back to the source
#it is critical not to use `_` in name of reference for compatibility with ddocent and freebayes

mv Tzo_scaffolds_TzC0402G_contam_R1R2_noIsolate.fasta ./reference.ssl.Tzo-C-0402G-R1R2-contam-noisolate.fasta


# 2. Mapping with the GenErode pipeline

Navigate to https://github.com/philippinespire/pire_taeniamia_zosterophora_lcwgs/tree/main/2nd_sequencing_run/GenErode_20k_2runs_2 to find GenErode config files and scripts.

Used GenErode v.0.6.0. https://github.com/NBISweden/GenErode

## GenErode for Tzo lcwgs

Ran pre-processing pipeline up to re-paired fastq files, will use these as input files for GenErode

## Set up your analysis folder

Make a copy of the template folder, renaming it according to your species name.

cp /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/scripts/GenErode_Wahab/GenErode_templatedir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2

Make directories within this directory to hold your config file, historical, modern, reference genome, and GERP scores.

mkdir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/config  
mkdir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/historical  
mkdir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/modern  
mkdir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/reference  
mkdir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/gerp_outgroups  
mkdir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/mitochondria  

All sequencing data and reference genomes have to be within the main analysis directory. Copy those to the appropriate subdirectories.

cp /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/fq_fp1_clmp_fp2_fqscrn_rprd/Tzo-A* /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/historical  
cp /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/fq_fp1_clmp_fp2_fqscrn_rprd/Tzo-C* /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/modern  
cp /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/1st_sequencing_run/fq_fp1_clmp_fp2_fqscrn_rprd/Tzo-A* /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/historical  
cp /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/1st_sequencing_run/fq_fp1_clmp_fp2_fqscrn_rprd/Tzo-C* /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/modern

Delete any .fq.gz files that have no reads left (re-pair output file would have listed these as errors). Files with no reads will cause GenErode to fail.

Subset reference genome to scaffold with 20K bases and above:  
README and script for removing small scaffolds is in the PSMC REU repo: https://github.com/philippinespire/REUs/tree/master/2022_REU/PSMC  
perl removesmalls.pl 20000 /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/refGenome/reference.ssl.Tzo-C-0402G-R1R2-contam-noisolate.fasta > reference.denovoSSL.Tzo20k.fasta

Check number of scaffolds of filtered assembly  
cat reference.denovoSSL.Tzo20k.fasta | grep "^>" | wc -l  
11,215 scaffolds left after filtering  
Non-filtered genome: 214,563 scaffolds

Check length of filtered assembly    
cat reference.denovoSSL.Tzo20k.fasta | grep -v "^>" | tr "\n" "\t" | sed 's/\t//g' | wc -c  
412,014,427  
Non-filtered genome: 935,413,606

Copy new reference genome  
cp /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/refGenome/reference.denovoSSL.Tzo20k.fasta /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/reference

Downloaded outgroup genomes from NCBI using wget (these need to be .gz files)

Leaving mitochondria directory empty for now as not planning on running mitochondrial analyses right away.

## Get Newick tree  

upload the list of species names to timetree.org ("Load a List of Species" at the bottom)  
download in Newick format

Upload newick tree file to your analysis folder in Wahab.

In the tree, rename the focal species with the name of the reference assembly file

## Edit the config files 

Copy the config file

home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/scripts/GenErode_Wahab/config /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/config 

Changes made to the config.yaml file:

ref_path: "/archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/reference/reference.denovoSSL.Tzo20k.fasta"

historical_samples: "config/tzo_historical_samples_2runs.txt" # leave empty ("") if not run for historical samples.  

modern_samples: "config/tzo_modern_samples_2runs.txt" # leave empty ("") if not run for modern samples. 

fastq_processing: False

bam_rmdup_realign_indels: True #This is needed to merge the sample bam files that were sequenced across multiple lanes

map_historical_to_mitogenomes: False

historical_bam_mapDamage: True

historical_rescaled_samplenames: Type out all historical sample names in this format: ["TzoABol008", "TzoABol009", "TzoABol010"]

Plink parameters:

Minimum SNP count. For example, 10, 25. Abbreviation in file name: homsnp.
homozyg-snp: 25  

Minimum size of ROH in kilobases. Abbreviation in file name: homkb.
homozyg-kb: 100  

Window size for ROH estimation. For example, 20, 50, 100, 1000. 
Abbreviation in file name: homwinsnp.
homozyg-window-snp: 250   

Maximum number of heterozygote sites per window. 
For example, 1 for a stringent analysis, 3 for a relaxed setting. 
Abbreviation in file name: homwinhet.
homozyg-window-het: 3  

Maximum number of missing sites per window. 
For example, 1 or 5 for a stringent analysis, 
10 for intermediate filtering and 15 for relaxed filtering. 
Abbreviation in file name: homwinmis.
homozyg-window-missing: 15  

Maximum number of heterozygote sites per ROH. 
For example, 1 for a stringent analysis, 
3 as a relaxed setting. Disable this parameter by setting it to 999. 
Abbreviation in file name: homhet.
homozyg-het: 750 

gerp: True

gerp_ref_path: "/archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/gerp_outgroups"

tree: "/archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/Tzo_outgroup_list.nwk"

min_gerp: 0  

max_gerp: 1000

## Prepare historical and modern config files
Prepared historical and modern config files in Google Sheets

File format + header: samplename_index_lane readgroup_id readgroup_platform path_to_R1_fastq_file path_to_R2_fastq_file. Note that the 3-underscore convention in sample_index_lane must be followed must remove hyphens/underscores from sample names. Make sure space-separated, not tab.

To get file from Google Sheets into the space-delimited .txt format for GenErode:  
Download from Google Sheets as a .tsv  
Open .tsv in Notepad ++  
Select all text, then click on “Search” then navigate to “Replace”  
Find what: \t  
Replace with: type in one space  
Search mode: Extended  
Click “Replace All”  
Save file as a “.txt” and transfer to Wahab  

## Run

Copy and run the sbatch file 

cp /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/scripts/GenErode_Wahab/run_GenErode.sbatch /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2

cd /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2

sbatch run_GenErode.sbatch 

Need GenErode to run until it generates realigned and rescaled bam files (modern files will only be realigned). These bam files will be used as the input for ANGSD.

# 3. Downstream analysis using ANGSD

Downstream analysis was run using ANGSD v. 0.940. 

Scripts are adapted from Nina Therkildsen's lab's genomic-data-analysis repository (https://github.com/therkildsen-lab/genomic-data-analysis).

Scripts are found in the Scripts directory and intermediate files in the Files directory unless otherwise noted.

## Start by identifying a list of SNPs for further analysis.

This SNP list was generated to ensure that analyses are run on the same loci across all 4 populations.

Unless otherwise noted, this list of SNPs was used for all subsequent analyses.

Filters used to generate a list of SNPs: minimum depth filters that corresponded to 1x per individual, maximum depth filters corresponding to 15x per individual, a minimum individual filter of 50% of the total number of individuals (111), a map quality filter of 30, a minimum allele frequency filter of 0.001, and a SNP p-value of 0.000001.

Transitions were excluded in all subsequent analyses.

sbatch snp_calling.sbatch <directory where your bam files are stored>

Yielded a list of ~1.78 million SNPs. 

## Get genotype likelihoods and generate a BEAGLE file for all populations.

Beagle file for all populations will be used for population structure (PCANGSD).

sbatch get_beagle.sbatch <directory where your bam files are stored>

Output: angsd_depth1_15_notrans.beagle.gz (found at 10.6084/m9.figshare.28832192 due to large file size)

Output run on neutral SNP list: angsd_allpop_neutral_alloutputs.beagle.gz (found at 10.6084/m9.figshare.28832192 due to large file size)

## PCANGSD- PCA and Admixture analyses 

Running PCANGSD for PCA and Admixture analysis. The MAP test yielded K=2 as most supported. Admixture was also set manually to run for K=3-5 (using the e= argument).

PCA: sbatch pcangsd_pca.sbatch <directory where your Beagle file is stored>

Output: ".cov" matrix is used to generate a PCA plot in the pca.R script.  
Output for all SNPs: angsd_notrans_snps_pca.cov  
Ouput for neutral SNPs: angsd_allpop_neutral_pca.cov  

Admixture: sbatch pcangsd_admix.sbatch <directory where your Beagle file is stored>

Output: ".Q" files for K=2-5 are used to generate admixture plots in the admixture.R script.  
Output for all SNPs and K=2: angsd_admix_notrans.admix.2.Q  
Output for neutral SNPs and K=2: angsd_allpop_neutral_default.admix.2.Q  

## Get genotype likelihoods, site allele frequencies, and minor allele frequencies for each population. 

Site allele frequencies (.saf) files are needed to calculate the SFS and genetic diversity metrics.

This script needs to be run for each population. A text file with a list of bam files for each population is needed in each script.
We utilized the reference genome for both -ref and -anc because an ancestral genome was not available. 

sbatch saf_beagle_maf.sbatch <directory with your bam files>

The .saf.idx files are used in the next step to generate a folded SFS for each population. They are also used in pairwise Fst calculations.

The saf outputs are: abol_sites_notrans.saf.idx, amta_sites_notrans.saf.idx, cbol_sites_notrans.saf.idx, cmta_sites_notrans.saf.idx

The saf outputs run on the neutral SNP list are: abol_neutral_notrans.saf.idx, amta_neutral_notrans.saf.idx, cbol_neutral_notrans.saf.idx, cmta_neutral_notrans.saf.idx

The .mafs.gz files are used in the selection.R script to run selection scans.

The maf outputs are: abol_sites_notrans.mafs.gz, amta_sites_notrans.mafs.gz, cbol_sites_notrans.mafs.gz, cmta_sites_notrans.mafs.gz 

The maf outputs run on the neutral SNP list are: abol_neutral_notrans.mafs.gz, amta_neutral_notrans.mafs.gz, cbol_neutral_notrans.mafs.gz, cmta_neutral_notrans.mafs.gz

## Generate the site frequency spectrum (SFS) for each population.
Generated a folded SFS because we did not have a known ancestral state genome.

This script needs to be run for each population. It utilizes the .saf.idx files generated in the last step.

sbatch sfs.sbatch <directory with your saf files>

Generates a .sfs file for each population that will be sused as an input for the saf2theta command in the next step.

## Calculate per-site thetas using the saf2theta command.

This script needs to be run for each population. It utilizes the .saf.idx and the .sfs files generated in the last steps.

sbatch saf2theta.sbatch <directory with your saf and sfs files>

Generates a thetas.idx file for each population that will be used in the next step.

## Calculate neutrality test statistics using the do_stat command.

This script needs to be run for each population. It utilizes the .thetas.idx file generated in the last step.

sbatch thetastat.sbatch <directory where your thetas.idx file is>

The output .thetas.idx.pestPG file is used for statistical analysis in the geneticdiversity.R script. Since we are using a folded SFS (unknown ancestral state), we are able to generate Watterson's theta (thetaW), nucleotide diversity (thetaD), and Tajima's D. 


## Calculate pairwise Fst 

This script needs to be run for each pair of populations. Each location/time period is one population (e.g. Malampaya Historical is one population). 
It requires the .saf.idx files for each population.

sbatch fst.sbatch <directory where the .saf.idx files are>

Unweighted and weighted Fst estimates are generated. We reported weighted estimates in our manuscript.


# 4. Statistical Analysis in R

## Identifying selection using the selection.R script 

Using the adapted CMH test to identify SNP candidates for selection at both sampling sites, and the adapted Chi-squared test to identify SNP candidates for selection 
at each sampling site individually. A generation time of one year was used. Effective density estimates come from the "Ne_estimation.R" script anduse Jorde and Ryman's method (2007).
SNPs with a false detection rate (FDR) less than 0.05 were considered to be under selection.

A list of all SNPs with the selection candidates excluded (i.e. putatively neutral SNPs) is compiled at the end of the selection.R script. ~1.3 million SNPs are on this neutral list.
This neutral SNP list is then used to run all downstream analyses in ANGSD to evaluate how selection may have impacted genetic diversity and population structure results. 
See the "Downstream analysis in ANGSD" section of the README to run these analyses. The following statistical analyses in R are run for both the 1.78 million SNP output files and the 1.3 million neutral SNP output files. 

Scripts are found in the Scripts directory and intermediate files in the Files directory unless otherwise noted.

## Generating PCA plots using the pca.R and pca_neutral.R scripts

Utilizing the .cov matrices from PCANGSD as the input. pca.R utilizes all 1.78 million SNPs, and pca_neutral.R utilizes the 1.3 million neutral SNPs. 


## Generating Admixture plots using the admixture.R and admixture_neutral.R scripts

Utilizing the ".Q" proportion outputs from PCANGSD as the input. admixture.R utilizes all 1.78 million SNPs, and admixture_neutral.R utilizes the 1.3 million neutral SNPs. 

## Estimating effective population size (Ne) using Ne_estimation.R and Ne_estimation_neutral.R

Using the .mafs.gz outputs from the "saf_beagle_maf.sbatch" ANGSD script. 

Code developed from Jorde & Ryman 2007 and the NeEstimator manual v.2.1

Ne_estimation.R is used to generate Ne estimates for the adapted CMH and adapted Chi-squared selection scan tests. 

Ne_estimation_neutral.R is used to generate Ne estimates from the 1.3 million neutral SNPs. These estimates are reported in our manuscript. 

## Analyzing changes in genetic diversity: Watterson's theta, nucleotide diversity (pi), and Tajima's D using the geneticdiversity.R and geneticdiversity_neutral.R scripts. 

Using the .thetas.idx.pestPG outputs from ANGSD.

Watterson's theta and nucleotide diversity were originally plotted against sequencing depth (mean depth per individual) to evaluate any depth based correlations that may be biasing results. 

This analysis identified that genetic diversity was sensitive to sequencing depth below 3x or above 6x (Figures S1-S2); we therefore restricted analyses on genetic diversity metrics to the 2,291 contigs with 3-6x depth. The following statistical analyses for all three metrics were run on this 3-6x depth range (452,496 SNPs). 




