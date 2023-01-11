#!/bin/bash
#SBATCH --job-name=Updating_GSA_NR03_allsamples
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --nodes=1

module load PLINK/1.9-beta6-20190617

/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/2020-Epithelial_Oslo_Deconvolution/ongoing/bed/GSA_update/update_build.sh \
/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03_allsamples/BED/CeDNN_NR03_allsamples \
/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/2020-Epithelial_Oslo_Deconvolution/ongoing/bed/GSA_script/GSAMD-24v1-0_20011747_A5-b37.strand \
/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03_allsamples/BED/GSA_converted/CeDNN_NR03_allsamples_converted
