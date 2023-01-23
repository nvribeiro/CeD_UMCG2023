#!/bin/bash
#SBATCH --job-name=merge_data_allsamples
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --nodes=1

module load R/4.2.2-foss-2022a-bare

Rscript merge_genotypes.R
