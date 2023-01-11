## Originally written by Ryan Schubert for the Wheeler Lab, github@RyanSchu Lab github@wheelerlab script to convert illumina final report 
## gentoype data to lgen data. This script will generate the .lgen file and .map file used for plink format. For instructions on how to 
## generate the .fam file please see the .md file

library(dplyr)
library(tidyr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--illumina", help="file path of the sample list")
parser$add_argument("-o", "--outputdir", help="output directory")
args <- parser$parse_args()

Finalreport<-as.data.frame(read.table(file=args$illumina, sep='\t', header = T))
Finalreport<-filter(Finalreport, Allele1...Forward != 'I')
Finalreport["empty"]="0"

map<-select(Finalreport, Chr, SNP.Name, empty, Position)
map<-map[!duplicated(map),]
map<-map[complete.cases(map),]

lgen<-select(Finalreport, Sample.Index, Sample.ID, SNP.Name, Allele1...Forward, Allele2...Forward)
lgen<-lgen[!duplicated(lgen),]
lgen<-lgen[complete.cases(lgen),]
lgen<-filter(lgen, Allele1...Forward != "-" & Allele2...Forward != "-")

write.table(map, file = paste(args$outputdir,"/CeDNN_NR03_allsamples.map",sep=""), sep = "\t", col.names = F, row.names = F, quote = F)
write.table(lgen, file = paste(args$outputdir,"/CeDNN_NR03_allsamples.lgen",sep=""), sep = "\t", col.names = F, row.names = F, quote = F)
