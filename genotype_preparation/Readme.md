# Readme
Created by:		 umcg-nribeiro (Nathan Vinicius Ribeiro)
Created on:		 Fri Jan 20 17:41:47 CET 2023
Contact at:		 n.v.ribeiro@umcg.nl
Principal investigator:  Iris H. Jonkers		

# Folders
scripts_genotypes

# Description
This folder contains the scripts necessary to convert a GSA Illumina FinalReport.txt to a .vcf ready to imputation using the Michigan Imputation Server.
It also has the scripts to check for quality control of the imputation and prepare the files to run demuxlet.

# Detailed information for each script
01_merge_genotypes.R	Merges the two CeDNN genotypes datasets
01_merge_genotypes.sh	Runs merge_genotypes.R
02_illumina_to_lgen.R	Converts a GSA Illumina FinalReport.txt to a pair of .lgen and .map files
02_illumina_to_lgen.sh	Runs illumina_to_lgen.R
03_creating_fam.R	Creates a .fam file to be used in combination with the .lgen and .map
03_create_fam.sh	Runs creating_fam.R
04_create_bed.sh	PLINK to create BED file using as input the .lgen, .map and .fam
05_update_GSA.sh	Update the the variants names from the GSA format to compatible with reference
06_make_freq.sh		Create a frequency file, needed to create the .vcf file
07_HRC_formatting.sh	Creates a .vcf file per chromosome in the HRC format taking a .bim file and .frq
08_check_vcf.sh		Check if the vcf files are ok to be sent to imputation
09_Compressing_VCFS.sh	Compress vcf files
10_Submitting_vcf_michigan.sh	Sends the vcf files for imputation in the Michigan server
11_post_imputation_qc.sh	Runs IC tool to check the quality of the imputation
12_prepare_vcf_demuxlet.sh	Filters the SNPs based on the QC, merge chromosomes in one vcf file, subsets for just the target samples (genotypes) and converts the chrmosomes ID to prepare for liftover
13_liftover_h19tohg38.sh	Converts the information from the hg19 genome to hg38 (used in the RNA data alignment)
14_subset_vcf_per_batch.sh	Creates one .vcf per batch containing the only the genotypes for that batch
15_liftover_per_batch.sh	Same as above but for my two batches separately and also filters just exons
16_demuxlet.sh			Runs demuxlet to demultiplex the samples based on genotype using the imputed vcf file
