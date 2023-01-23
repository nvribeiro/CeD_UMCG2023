#!/bin/bash
#SBATCH --job-name=demuxlet_A2_160123
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16gb
#SBATCH --nodes=1

INPUT='/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/processed/count/NR03_A2_count_161222_4/outs/possorted_genome_bam.bam'
VCF='/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/imputed_unzipped/vcf_filtered/NR03_imputed_liftedover_onlyexons.vcf'
OUTPUT='/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/demuxlet'

module load demuxlet/20200417-3ab507c-GCCcore-7.3.0

demuxlet --sam $INPUT --vcf $VCF --field GT --out ${OUTPUT}/A2_demuxlet
