#!/bin/bash

#SBATCH --job-name=NR03_fam_file_all_samples
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --nodes=1

module load R/4.2.2-foss-2022a-bare
module list

Rscript creating_fam.R \
--input /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03_allsamples/LGEN/CeDNN_NR03_allsamples.lgen \
--outputdir /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03_allsamples/LGEN
