#!/bin/bash
#SBATCH --job-name=subset_vcf_per_batch
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --nodes=1

module load BCFtools/1.16-GCCcore-11.3.0

# Subseting the vcf just with the samples I want
INPUT='/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/processed/genotypes/CeDNN_NR03_imputation/imputed_unzipped'
SAMPLES_A=30_219-0607,51_219-1056,13_CeDNN-A0325,20_CeDNN-A0339,21_CeDNN-A0340,23_CeDNN-A0357,24_CeDNN-A0381
SAMPLES_F=51_219-1056,13_CeDNN-A0325,15_CeDNN-A0326,19_CeDNN-A0338,21_CeDNN-A0340

bcftools view -Ov -s ${SAMPLES_A} ${INPUT}/vcf_filtered/NR03_imputed_merged.vcf.gz > ${INPUT}/vcf_filtered/NR03_imputed_mysamples_A.vcf
bcftools view -Ov -s ${SAMPLES_F} ${INPUT}/vcf_filtered/NR03_imputed_merged.vcf.gz > ${INPUT}/vcf_filtered/NR03_imputed_mysamples_F.vcf

# Converting the chromosome IDs to chr
bcftools annotate --rename-chrs ${INPUT}/vcf_filtered/chr_name_conv.txt ${INPUT}/vcf_filtered/NR03_imputed_mysamples_A.vcf -Ov -o ${INPUT}/vcf_filtered/NR03_imputed_mysamples_A_with_chr.vcf
bcftools annotate --rename-chrs ${INPUT}/vcf_filtered/chr_name_conv.txt ${INPUT}/vcf_filtered/NR03_imputed_mysamples_F.vcf -Ov -o ${INPUT}/vcf_filtered/NR03_imputed_mysamples_F_with_chr.vcf
