#!/bin/bash
#SBATCH --job-name=post_imputation_QC_NR03
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=1:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --nodes=1

module load PerlPlus/5.30.0-GCCcore-7.3.0-v19.08.1
module load libgd/2.2.5-GCCcore-7.3.0-lor

export PERL5LIB=~/lib/perl5/site_perl:$PERL5LIB

perl /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/tools/IC/ic.pl \
-d /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/imputed_unzipped \
-r /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/tools/HRC.r1-1.GRCh37.wgs.mac5.sites.tab \
-o /groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/post_imputation_QC
