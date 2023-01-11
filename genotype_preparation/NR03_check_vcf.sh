#!/bin/bash
#SBATCH --job-name=NR03_check_vcf_allsamples
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --nodes=1

cd /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03_allsamples/checkVCF

for chr in {1..23}; do \
python /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/tools/checkVCF.py \
-r /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/tools/hs37d5.fa \
-o out-chr${chr} /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03_allsamples/HRC_formated/CeDNN_NR03_allsamples_converted-updated-chr${chr}.vcf;
done
