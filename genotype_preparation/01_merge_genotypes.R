setwd('/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/genotypes/raw/all_samples')
library(dplyr)

CeD20 <- as.data.frame(read.table('CeDNN_2020_FinalReport.txt', sep='\t', skip = 9, header = T))
CeD19 <- as.data.frame(read.table('CeDNN_GSAMD1_21112019_FinalReport.txt', sep='\t', skip = 9, header = T))

new_dataset <- bind_rows(CeD20, CeD19)
write.table(new_dataset, file = 'CeD_NR03_allsamples.txt', sep = "\t", col.names = T, row.names = F, quote = F)
