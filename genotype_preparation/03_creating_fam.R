# Creating .fam file for PLINK using .lgen file as input
library(dplyr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--input", help="file path of the sample list")
parser$add_argument("-o", "--outputdir", help="output directory")
args <- parser$parse_args()

lgen <- as.data.frame(read.table(file = args$input,
                   header = FALSE,
                   sep = "\t"))

names(lgen) <- c('FID', 'IID', 'SNP', 'A1', 'A2')

fam <- data.frame(FID = lgen$FID, IID = lgen$IID, Paternal_ID = 0,
                  Maternal_ID = 0, Sex = 0, Phenotype = 0)
fam <- fam %>% distinct(IID, .keep_all = TRUE)

write.table(fam, file = paste0(args$outputdir,"/CeD.fam"), sep = "\t", col.names = F, row.names = F, quote = F)
