#!/bin/bash
#SBATCH --job-name=NR03_make_bed_050123
#SBATCH --output=../%x.out
#SBATCH --error=../%x.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --nodes=1

module load PLINK/1.9-beta6-20190617

plink --lfile /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03/LGEN/CeDNN_NR03 \
--make-bed --out /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03/BED/CeDNN_NR03
