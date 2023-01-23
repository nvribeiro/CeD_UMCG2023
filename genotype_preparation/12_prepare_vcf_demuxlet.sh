#!/bin/bash
#SBATCH --job-name=filter_vcf
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --nodes=1

INPUT='/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/processed/genotypes/CeDNN_NR03_imputation/imputed_unzipped'

module load BCFtools/1.16-GCCcore-11.3.0

# Excluding SNPs with R2 < 0.8 and MAF < 0.02
for chr in {1..23}; do \
bcftools filter -Oz -e 'INFO/R2<0.8 && INFO/MAF<0.02' ${INPUT}/chr${chr}.dose.vcf.gz > ${INPUT}/vcf_filtered/chr${chr}_filtered.vcf.gz ; \
done

# Merging all chromosomes into one vcf file
bcftools concat -f ${INPUT}/vcf_filtered/sample_list.txt -Oz -o ${INPUT}/vcf_filtered/NR03_imputed_merged.vcf.gz

# Subseting the vcf just with the samples I want
SAMPLES=30_219-0607,51_219-1056,13_CeDNN-A0325,15_CeDNN-A0326,19_CeDNN-A0338,20_CeDNN-A0339,23_CeDNN-A0357,24_CeDNN-A0381
bcftools view -Ov -s ${SAMPLES} ${INPUT}/vcf_filtered/NR03_imputed_merged.vcf.gz > ${INPUT}/vcf_filtered/NR03_imputed_mysamples.vcf

# Converting the chromosome IDs to chr
bcftools annotate --rename-chrs ${INPUT}/vcf_filtered/chr_name_conv.txt ${INPUT}/vcf_filtered/NR03_imputed_mysamples.vcf -Ov -o ${INPUT}/vcf_filtered/NR03_imputed_mysamples_with_chr.vcf
