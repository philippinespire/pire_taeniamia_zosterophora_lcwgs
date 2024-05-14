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

##Subset reference genome to scaffold with 20K bases and above:
README and script for removing small scaffolds is in the PSMC REU repo: https://github.com/philippinespire/REUs/tree/master/2022_REU/PSMC
perl removesmalls.pl 20000 /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/refGenome/reference.ssl.Tzo-C-0402G-R1R2-contam-noisolate.fasta > reference.denovoSSL.Tzo20k.fasta

#Check number of scaffolds of filtered assembly
cat reference.denovoSSL.Tzo20k.fasta | grep "^>" | wc -l
11,215 scaffolds left after filtering
Non-filtered genome: 214,563 scaffolds

#Check length of filtered assembly:
cat reference.denovoSSL.Tzo20k.fasta | grep -v "^>" | tr "\n" "\t" | sed 's/\t//g' | wc -c
412,014,427
Non-filtered genome: 935,413,606

#Copy new reference genome
cp /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/refGenome/reference.denovoSSL.Tzo20k.fasta /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/reference

Downloaded outgroup genomes from NCBI using wget (these need to be .gz files)

Leaving mitochondria directory empty for now as not planning on running mitochondrial analyses right away.

##Get Newick tree

Get Newick tree:

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


##Prepare historical and modern config files
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

##Need GenErode to run until it generates a .ancestral.rates.gz file (should be under results/gerp/referencename/)
