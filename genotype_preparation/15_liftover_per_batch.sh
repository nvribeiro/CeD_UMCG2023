#!/bin/bash
#SBATCH --job-name=liftover_NR03_per_batch
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --nodes=1

FOLDER=/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/processed/genotypes/CeDNN_NR03_imputation/imputed_unzipped/vcf_filtered

module load picard/2.26.10-Java-8-LTS
module load BEDTools/2.30.0-GCCcore-11.3.0
module load BCFtools/1.16-GCCcore-11.3.0

# A batch
java -jar ${EBROOTPICARD}/picard.jar LiftoverVcf \
I=${FOLDER}/NR03_imputed_mysamples_A_with_chr.vcf \
O=${FOLDER}/NR03_imputed_liftedover_A.vcf \
CHAIN=/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/tools/hg19ToHg38.over.chain \
REJECT=${FOLDER}/NR03_rejected_variants.vcf \
R=/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/GRCh38_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
WARN_ON_MISSING_CONTIG=true \
RECOVER_SWAPPED_REF_ALT=true

# F batch
java -jar ${EBROOTPICARD}/picard.jar LiftoverVcf \
I=${FOLDER}/NR03_imputed_mysamples_F_with_chr.vcf \
O=${FOLDER}/NR03_imputed_liftedover_F.vcf \
CHAIN=/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/tools/hg19ToHg38.over.chain \
REJECT=${FOLDER}/NR03_rejected_variants.vcf \
R=/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/GRCh38_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
WARN_ON_MISSING_CONTIG=true \
RECOVER_SWAPPED_REF_ALT=true

# Filter INFO > 0.9
bcftools filter -Oz -e 'INFO/R2<0.9' ${FOLDER}/NR03_imputed_liftedover_A.vcf > ${FOLDER}/NR03_imputed_liftedover_A_INFO09.vcf
bcftools filter -Oz -e 'INFO/R2<0.9' ${FOLDER}/NR03_imputed_liftedover_F.vcf > ${FOLDER}/NR03_imputed_liftedover_F_INFO09.vcf

# Keep only exons
bedtools intersect -a ${FOLDER}/NR03_imputed_liftedover_A_INFO09.vcf -b /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/tools/hg38_exons.bed -wa > ${FOLDER}/NR03_imputed_liftedover_A_onlyexons.vcf
bedtools intersect -a ${FOLDER}/NR03_imputed_liftedover_F_INFO09.vcf -b /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/tools/hg38_exons.bed -wa > ${FOLDER}/NR03_imputed_liftedover_F_onlyexons.vcf
