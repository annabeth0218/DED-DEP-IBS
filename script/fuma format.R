library(data.table)
library(readxl)
library(BSgenome.Hsapiens.UCSC.hg38) # using GRCh38/hg38
library(tidyverse)

tp <- gwas.ibs |>
  filter(str_detect(SNP, "^rs")) |>
  select(SNP, CHR, BP, A1, A2, NMISS, OR, SE, L95, U95, STAT, P) 
write.table(tp, "gwas_ibs.txt", 
            sep = " ", quote = FALSE, row.names = FALSE)

# top
top_200 <- data.table(mtag.ibs)
top_200 <- top_200[order(mtag_pval)][1:200]
write_csv(top_200, 'mtag_gwas_snp/one-mtag-ibs.csv')


# pos + eQTL + ci = genes.disease.1g/1m
f.ibs.1m <- read_excel("fuma/snp2gene.xlsx", sheet = "ibs-1m")
f.ibs.1g <- read_excel("fuma/snp2gene.xlsx", sheet = "ibs-1g")
f.ded.1m <- read_excel("fuma/snp2gene.xlsx", sheet = "ded-1m")
f.ded.1g <- read_excel("fuma/snp2gene.xlsx", sheet = "ded-1g")
f.dep.1m <- read_excel("fuma/snp2gene.xlsx", sheet = "dep-1m")
f.dep.1g <- read_excel("fuma/snp2gene.xlsx", sheet = "dep-1g")
