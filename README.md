# pire_taeniamia_zosterophora_lcwgs
Repository for bioinformatics and analysis of historical and contemporary low-coverage whole genome sequencing data of Taeniamia zosterophora collected in the Philippines.

# Pre-processing
Following the pire_fq_gz_processing repo

Files were already downloaded, and the decode file looks fine, so moving to Step 5: renaming

# 5. Perform a renaming dry-run
   Run renameFQGZ.bash to view the original and new file names and create tsv files to store the original and new file naming conventions.
 
log into a compute node interactively so this goes faster
salloc

once you have the compute node, procede
bash /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/fq_raw/renameFQGZ.bash Tzo-LCWGS-FullSeq-SequenceNameDecode.tsv 

Dry run looks okay, so will proceed to renaming for real.

# 6. Rename the files for real

bash /archive/carpenterlab/pire/pire_taeniamia_zosterophora_lcwgs/2nd_sequencing_run/fq_raw/renameFQGZ.bash Tzo-LCWGS-FullSeq-SequenceNameDecode.tsv rename

#you will need to say y 2X
