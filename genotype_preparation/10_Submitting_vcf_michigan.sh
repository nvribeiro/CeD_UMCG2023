#!/bin/bash
#SBATCH --job-name=Submitting_VCFs
#SBATCH --output=../script_outputs/%x.out
#SBATCH --error=../script_outputs/%x.err
#SBATCH --time=1:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --nodes=1

module load cURL/7.83.0-GCCcore-11.3.0;

TOKEN="eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJuYW1lIjoiQWFyb24gUmFtaXJleiIsImFwaSI6dHJ1ZSwibWFpbCI6ImFhcm9uZHJzMkBob3RtYWlsLmNvbSIsImV4cGlyZSI6MTY3NTg0Njk0NzQ0NCwidXNlcm5hbWUiOiJhYXJvbmRyczIifQ.7GNRu220fs4OrHWWKqslJbW2fPx9G8RPmuT7-2UabAU"

folder="/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/converted/CeDNN_NR03/vcf_compressed_for_imputation"

curl https://imputationserver.sph.umich.edu/api/v2/jobs/submit/minimac4 \
  -H "X-Auth-Token: ${TOKEN}" \
  -F "files=@${folder}/chr1.vcf.gz" \
  -F "files=@${folder}/chr2.vcf.gz" \
  -F "files=@${folder}/chr3.vcf.gz" \
  -F "files=@${folder}/chr4.vcf.gz" \
  -F "files=@${folder}/chr5.vcf.gz" \
  -F "files=@${folder}/chr6.vcf.gz" \
  -F "files=@${folder}/chr7.vcf.gz" \
  -F "files=@${folder}/chr8.vcf.gz" \
  -F "files=@${folder}/chr9.vcf.gz" \
  -F "files=@${folder}/chr10.vcf.gz" \
  -F "files=@${folder}/chr11.vcf.gz" \
  -F "files=@${folder}/chr12.vcf.gz" \
  -F "files=@${folder}/chr13.vcf.gz" \
  -F "files=@${folder}/chr14.vcf.gz" \
  -F "files=@${folder}/chr15.vcf.gz" \
  -F "files=@${folder}/chr16.vcf.gz" \
  -F "files=@${folder}/chr17.vcf.gz" \
  -F "files=@${folder}/chr18.vcf.gz" \
  -F "files=@${folder}/chr19.vcf.gz" \
  -F "files=@${folder}/chr20.vcf.gz" \
  -F "files=@${folder}/chr21.vcf.gz" \
  -F "files=@${folder}/chr22.vcf.gz" \
  -F "files=@${folder}/chr23.vcf.gz" \
  -F "refpanel=hrc-r1.1" \
  -F "population=eur"
