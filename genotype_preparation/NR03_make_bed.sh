#!/bin/bash
#SBATCH --job-name=NR03_make_bed_allsamples
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=82gb
#SBATCH --nodes=1

module load PLINK/1.9-beta6-20190617

plink --lfile /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03_allsamples/LGEN/CeDNN_NR03_allsamples \
--make-bed --out /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03_allsamples/BED/CeDNN_NR03_allsamples
