#!/bin/bash
#SBATCH --job-name=NR03_HRC_formatting_allsamples
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=1:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb
#SBATCH --nodes=1

module load PerlPlus/5.30.0-GCCcore-7.3.0-v19.08.1
module load PLINK/1.9-beta6-20190617

cd /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03_allsamples/BED/GSA_converted

perl /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/2020-Epithelial_Oslo_Deconvolution/ongoing/HRC-1000G-check-bim-v4/HRC-1000G-check-bim-NoReadKey.pl \
-b CeDNN_NR03_allsamples_converted.bim \
-f CeDNN_NR03_allsamples_converted.frq \
-r /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/2020-Epithelial_Oslo_Deconvolution/ongoing/HRC-1000G-check-bim-v4/tab/HRC.r1-1.GRCh37.wgs.mac5.sites.tab \
-h

sh Run-plink.sh
