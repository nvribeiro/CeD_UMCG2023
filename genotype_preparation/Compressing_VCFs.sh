#!/bin/bash
#SBATCH --job-name=Compressing_VCFs_allsamples
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=1:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --nodes=1

module load BCFtools

for chr in {1..23}; do \
bcftools view -O z -o "../vcf_compressed_for_imputation/chr${chr}.vcf.gz" "../vcf_for_imputation/chr${chr}.vcf"; \
done
