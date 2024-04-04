## GenErode for Tzo lcwgs

Ran pre-processing pipeline up to re-paired fastq files, will use these as input files for GenErode

## Set up your analysis folder

Make a copy of the template folder, renaming it according to your species name.

cp /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/scripts/GenErode_Wahab/GenErode_templatedir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode


Make directories within this directory to hold your config file, historical, modern, reference genome, and GERP scores.

mkdir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode/config
mkdir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode/historical
mkdir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode/modern
mkdir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode/reference
mkdir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode/gerp_outgroups
mkdir /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode/mitochondria

All sequencing data and reference genomes have to be within the main analysis directory. Copy those to the appropriate subdirectories.

cp /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/fq_fp1_clmp_fp2_fqscrn_rprd/Tzo-A* /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode/historical
cp /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/fq_fp1_clmp_fp2_fqscrn_rprd/Tzo-C* /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode/modern
cp /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/refGenome/reference.ssl.Tzo-C-0402G-R1R2-contam-noisolate.fasta /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode/reference

Downloaded outgroup genomes from NCBI using wget

Leaving mitochondria directory empty for now as not planning on running mitochondrial analyses right away.

## Edit the config files 

Copy the config file

home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/scripts/GenErode_Wahab/config /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode/config 

Changes made to the config.yaml file:

ref_path: "/archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode/reference/reference.ssl.Tzo-C-0402G-R1R2-contam-noisolate.fasta"

historical_samples: "config/tzo_historical_samples.txt" # leave empty ("") if not run for historical samples.
modern_samples: "config/tzo_modern_samples.txt" # leave empty ("") if not run for modern samples. 

fastq_processing: False

bam_rmdup_realign_indels: True #This is needed to merge the sample bam files that were sequenced across multiple lanes

map_historical_to_mitogenomes: False

historical_bam_mapDamage: True

historical_rescaled_samplenames: ["TzoABol008", "TzoABol009", "TzoABol010", "TzoABol020", "TzoABol023", "TzoABol026", "TzoABol030", "TzoABol034", "TzoABol036", "TzoABol038", "TzoABol043", "TzoABol044", "TzoABol045", "TzoABol047", "TzoABol050", "TzoABol052", "TzoABol053", "TzoABol054", "TzoABol055", "TzoABol056", "TzoABol057", "TzoABol058", "TzoABol059"
"TzoABol060", "TzoABol061", "TzoABol062", "TzoABol065", "TzoABol066", "TzoABol067", "TzoABol068", "TzoABol069", "TzoABol072", "TzoABol073", "TzoABol075", "TzoABol076", "TzoABol077", "TzoABol078", "TzoABol079", "TzoABol080", "TzoABol081", "TzoABol082", "TzoABol083", "TzoABol084", "TzoABol085", "TzoABol086", "TzoABol087", "TzoABol091", "TzoABol092", "TzoAMta001",
"TzoAMta002", "TzoAMta003", "TzoAMta004", "TzoAMta005", "TzoAMta006", "TzoAMta007", "TzoAMta008", "TzoAMta009", "TzoAMta010", "TzoAMta011", "TzoAMta012", "TzoAMta013", "TzoAMta014", "TzoAMta015", "TzoAMta016", "TzoAMta017", "TzoAMta018", "TzoAMta019", "TzoAMta020", "TzoAMta021", "TzoAMta022", "TzoAMta023", "TzoAMta024", "TzoAMta025", "TzoAMta026", "TzoAMta027",
"TzoAMta028", "TzoAMta029", "TzoAMta030", "TzoAMta031", "TzoAMvi001","TzoAMvi002", "TzoAMvi003", "TzoAMvi004", "TzoAMvi005", "TzoAMvi006", "TzoAMvi007", "TzoAMvi008", "TzoAMvi009", "TzoAMvi010", "TzoAMvi011", "TzoAMvi012", "TzoAMvi013", "TzoAMvi014", "TzoAMvi015", "TzoAMvi016", "TzoAMvi017", "TzoAMvi018", "TzoAMvi019", "TzoAMvi020", "TzoAMvi021", "TzoAMvi022",
"TzoAMvi023", "TzoAMvi024", "TzoAMvi025", "TzoAMvi026", "TzoAMvi027", "TzoAMvi028", "TzoAMvi029", "TzoAMvi030", "TzoAMvi031", "TzoAMvi032", "TzoAMvi033", "TzoAMvi034", "TzoAMvi035", "TzoAMvi036", "TzoAMvi037", "TzoAMvi038", "TzoAMvi039", "TzoAMvi040", "TzoAMvi041", "TzoAMvi042", "TzoAMvi043", "TzoAMvi044", "TzoAMvi045", "TzoAMvi046", "TzoAMvi047", "TzoAMvi048"]

Plink parameters:

# Minimum SNP count. For example, 10, 25. Abbreviation in file name: homsnp.
homozyg-snp: 25
# Minimum size of ROH in kilobases. Abbreviation in file name: homkb.
homozyg-kb: 100 
# Window size for ROH estimation. For example, 20, 50, 100, 1000. 
# Abbreviation in file name: homwinsnp.
homozyg-window-snp: 250 
# Maximum number of heterozygote sites per window. 
# For example, 1 for a stringent analysis, 3 for a relaxed setting. 
# Abbreviation in file name: homwinhet.
homozyg-window-het: 3 
# Maximum number of missing sites per window. 
# For example, 1 or 5 for a stringent analysis, 
# 10 for intermediate filtering and 15 for relaxed filtering. 
# Abbreviation in file name: homwinmis.
homozyg-window-missing: 15 
# Maximum number of heterozygote sites per ROH. 
# For example, 1 for a stringent analysis, 
# 3 as a relaxed setting. Disable this parameter by setting it to 999. 
# Abbreviation in file name: homhet.
homozyg-het: 750 

gerp: True

gerp_ref_path: "/archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode/gerp_outgroups"

tree: "/archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode/Tzo_outgroup_list.nwk"

min_gerp: 0
max_gerp: 1000


Prepared historical and modern config files in Google Sheets

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

cp /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/scripts/GenErode_Wahab/run_GenErode.sbatch /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode

cd /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode

sbatch run_GenErode.sbatch 
