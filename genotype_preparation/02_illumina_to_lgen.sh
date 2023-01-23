#!/bin/bash
#SBATCH --job-name=NR03_illumina_to_lgen
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --nodes=1

module load R/4.2.2-foss-2022a-bare

Rscript illumina_to_lgen.R \
--illumina  /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/raw/filtered_samples/CeD_NR03.txt \
--outputdir /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03
