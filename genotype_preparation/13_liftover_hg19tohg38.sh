#!/bin/bash
#SBATCH --job-name=liftover_NR03
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --nodes=1

FOLDER=/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/imputed_unzipped/vcf_filtered

module load picard/2.26.10-Java-8-LTS

# Liftover genome hg19 to hg38
java -jar ${EBROOTPICARD}/picard.jar LiftoverVcf \
I=${FOLDER}/NR03_imputed_mysamples_with_chr.vcf \
O=${FOLDER}/NR03_imputed_liftedover.vcf \
CHAIN=/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/tools/hg19ToHg38.over.chain \
REJECT=${FOLDER}/NR03_rejected_variants.vcf \
R=/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/GRCh38_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
WARN_ON_MISSING_CONTIG=true \
RECOVER_SWAPPED_REF_ALT=true
