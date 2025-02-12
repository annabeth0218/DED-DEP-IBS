library(data.table)
library(readxl)
library(BSgenome.Hsapiens.UCSC.hg38) # using GRCh38/hg38
library(tidyverse)

# import
gwas.dep <- read.csv('table/gwas_snp/with_A2/one_sample_gwas/depression_log.csv')
gwas.ded <- read.csv('table/gwas_snp/with_A2/one_sample_gwas/dryeye_case60_log.csv')
gwas.ibs <- read.csv('table/gwas_snp/with_A2/one_sample_gwas/ibs_case50_log.csv')

tp <- gwas.dep |>
  filter(str_detect(SNP, "^rs")) |>
  select(CHR, SNP, BP, A1, A2, NMISS, P) 
write.table(tp, "fuma_dep.txt", 
            sep = " ", quote = FALSE, row.names = FALSE)

tp <- overlap |>
  filter(str_detect(SNP, "^rs")) |>
  select(SNP, A1, NMISS, OR, SE, P) 
write.table(tp, "overlap_ded-dep.txt", 
            sep = " ", quote = FALSE, row.names = FALSE)

# top
top_200 <- data.table(mtag.ibs)
top_200 <- top_200[order(mtag_pval)][1:200]
write_csv(top_200, 'mtag_gwas_snp/one-mtag-ibs.csv')

# ANNOVAR res
tp <- read.table('fuma/dep-one-m_524702/FUMA_job524702/snps.txt', header = TRUE)
dep.1m <- tp |> filter(!is.na(gwasP))
tp <- dep.1m |>
  group_by(nearestGene) |>
  filter(gwasP == min(gwasP, na.rm = TRUE)) |>
  ungroup()
write.table(tp, 'fuma/dep-1m.txt',
            sep = "\t", quote = FALSE, row.names = FALSE)

# pos + eQTL + ci = genes.disease.1g/1m
f.ibs.1m <- read_excel("fuma/snp2gene.xlsx", sheet = "ibs-1m")
f.ibs.1g <- read_excel("fuma/snp2gene.xlsx", sheet = "ibs-1g")
f.ded.1m <- read_excel("fuma/snp2gene.xlsx", sheet = "ded-1m")
f.ded.1g <- read_excel("fuma/snp2gene.xlsx", sheet = "ded-1g")
f.dep.1m <- read_excel("fuma/snp2gene.xlsx", sheet = "dep-1m")
f.dep.1g <- read_excel("fuma/snp2gene.xlsx", sheet = "dep-1g")

ibs.dep.1m <- merge(f.ibs.1m, f.dep.1m, by.x = "ensg", by.y = "ensg")

tp <- f.ibs.1m |> 
  filter(!is.na(f.ibs.1m$pLI) & pLI != "NA") |>
  mutate(pLI_numeric = as.numeric(pLI)) |>
  slice_max(pLI_numeric, n = 20)
View(tp)

tp <- f.ibs.1m |> 
  filter(!is.na(f.ibs.1m$ncRVIS) & ncRVIS != "NA") |>
  mutate(ncRVIS_numeric = as.numeric(ncRVIS)) |>
  slice_max(ncRVIS_numeric, n = 20)
View(tp)