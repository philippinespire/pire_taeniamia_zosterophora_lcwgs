## ATLAS

## Recalibration

# move to your species directory - this directory should also contain your GenErode analysis directory

cd /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/

# make an ATLAS analysis directory

mkdir ATLAS_nurecal

# get sites with scores in each particular bin

gunzip -c GenErode_20k_2runs_2/results/gerp/reference.denovoSSL.Tzo20k.ancestral.rates.gz | awk -v OFS='\t' '($4 >= 1.5) {print $1,$2,$2}' > ATLAS_nurecal/gerp15.bed

gunzip -c GenErode_20k_2runs_2/results/gerp/reference.denovoSSL.Tzo20k.ancestral.rates.gz | awk -v OFS='\t' '($4 >= 1 && $4 < 1.5) {print $1,$2,$2}' > ATLAS_nurecal/gerp1.bed

gunzip -c GenErode_20k_2runs_2/results/gerp/reference.denovoSSL.Tzo20k.ancestral.rates.gz | awk -v OFS='\t' '($4 >= 0.5 && $4 < 1) {print $1,$2,$2}' > ATLAS_nurecal/gerp5.bed

# enter ATLAS directory

cd ATLAS_nurecal

# this is a good time to check the number of sites in your gerp*.bed files (`wc -l gerp15.bed` etc)... ATLAS wants something on the order of ~5 million sites total for its recalibration procedure.

Results: 96,410 gerp15.bed; 1,551,798 gerp1.bed; 7,034,373 gerp5.bed

# merge adjacent sites into regions using bedtools

module load bedtools

crun bedtools merge -i gerp15.bed > gerp15.merge.bed 

crun bedtools merge -i gerp1.bed > gerp1.merge.bed 

crun bedtools merge -i gerp5.bed > gerp5.merge.bed 

# remove singleton sites from merged .bed files.

awk -v OFS='\t' '!($2==$3)' < gerp15.merge.bed > gerp15.merge.nosingle.bed 

awk -v OFS='\t' '!($2==$3)' < gerp1.merge.bed > gerp1.merge.nosingle.bed 

awk -v OFS='\t' '!($2==$3)' < gerp5.merge.bed > gerp5.merge.nosingle.bed


#Get scripts
cp /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/scripts/ATLAS_wahab/atlas_recal_readuntilbeds_array.sbatch .
cp /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/scripts/ATLAS_wahab/atlas_recal_readuntilbeds_array.bash .


Modify the atlas_recal .bash script to reflect the path for the .sbatch script.
Modify the .sbatch script to reflect the number of .bed files you want to use and their names. You can also change the maximum depth - currently set to only use 10 reads nano maximum per site, which should be sufficient to calculate error.

#Create output directory
mkdir atlas_recal_depth10_update


#Copy bed files into output directory
Cp *.bed atlas_recal_depth10_update


#Run
bash atlas_recal_readuntilbeds_array.bash [directory with .bam files] [reference genome location] [output directory]


bash atlas_recal_readuntilbeds_array.bash /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/mkBAM  /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/GenErode_20k_2runs_2/reference/reference.denovoSSL.Tzo20k.fasta /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/ATLAS_nurecal/atlas_recal_depth10_update


## Calculate theta (heterozygosity)

#Create output directory

mkdir theta_recalibrated

Modify the atlas_recal .bash script to reflect the path for the .sbatch script.

#Calculating theta across the whole-genome using the --genomeWide command 

bash atlas_theta_albrecal_array.bash [directory with recalibrated output files] [output directory]

bash atlas_theta_albrecal_array.bash /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/ATLAS_nurecal/atlas_recal_depth10_update /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/ATLAS_nurecal/theta_recal_update


#Generate theta box-plots

#Download the *theta.txt.gz files to your personal computer.

#Then run the wrangle_plot_theta_template.R script (found in: pire_lcwgs_data_processing/scripts/ATLAS_wahab) locally in R. 


##Generating genotype likelihoods

#Creates a .glf.gz file for each bam file

mkdir glf

bash atlas_glf_albrecal_array.bash [directory with recalibrated bam and .json files] [output directory]
bash atlas_glf_albrecal_array.bash /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/ATLAS_nurecal/atlas_recal_depth10_update /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/ATLAS_nurecal/glf


##Generating MajorMinor VCF file:

#Creates one vcf file

#Creating text file with the names of all glf files for the samples between 0.1-5x depth: glf_0.1_5.txt, store this file in the same directory as .glf.gz files

#Using default method of MLE, minMAF=0.001 (from recommendation of setting minMAF lower than 1/2*# of individuals, so 1/199*2=0.0025, setting lower at 0.001)

sbatch atlas_majorminor_albrecal.sbatch [directory where .glf.gz files are stored]

sbatch atlas_majorminor_albrecal.sbatch /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/ATLAS_nurecal/glf 


##Converting VCF to Beagle file for use in ANGSD:

sbatch atlas_convertvcf_albrecal.sbatch [directory where .vcf.gz file is stored]

sbatch atlas_convertvcf_albrecal.sbatch /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/ATLAS_nurecal/glf 

#Now have one Beagle file, will copy this file to angsd directory (/archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/angsd) and start a 
#README for analysis in ANGSD and PCANGSD in this directory

