#!/bin/bash
#SBATCH --job-name=NR03_makefreq_allsamples
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --nodes=1

module load PLINK/1.9-beta6-20190617

plink --freq \
--bfile /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03_allsamples/BED/GSA_converted/CeDNN_NR03_allsamples_converted \
--out /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03_allsamples/BED/GSA_converted/CeDNN_NR03_allsamples_converted
