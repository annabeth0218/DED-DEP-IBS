library(data.table)
library(readxl)
library(BSgenome.Hsapiens.UCSC.hg38) # using GRCh38/hg38
library(tidyverse)

# import
gwas.dep <- read.csv('mtag_gwas_snp/onesample/gwas/depression_log_ADD.csv')
gwas.ded <- read.csv('mtag_gwas_snp/onesample/gwas/dryeye_case60_log_ADD.csv')
gwas.ibs <- read.csv('mtag_gwas_snp/onesample/gwas/ibs_case50_log_ADD.csv')

mtag.dep <- read.csv('mtag_gwas_snp/onesample/mtag/one_mtag_dep.csv')
mtag.ded <- read.csv('mtag_gwas_snp/onesample/mtag/one_mtag_dry60.csv')
mtag.ibs <- read.csv('mtag_gwas_snp/onesample/mtag/one_mtag_ibs50.csv')

gwas.dep.two <- read.csv('mtag_gwas_snp/twosample/gwas/two_dep_log_ADD.csv')
gwas.ded.two <- read.csv('mtag_gwas_snp/twosample/gwas/two_dry_case60_log_ADD.csv')
gwas.ibs.two <- read.csv('mtag_gwas_snp/twosample/gwas/two_ibs_case50_log_ADD.csv')

mtag.dep.two <- read.csv('mtag_gwas_snp/twosample/mtag/two_mtag_dep.csv')
mtag.ded.two <- read.csv('mtag_gwas_snp/twosample/mtag/two_mtag_dry60.csv')
mtag.ibs.two <- read.csv('mtag_gwas_snp/twosample/mtag/two_mtag_ibs50.csv')

format_snp_osmr <- read_excel("format-snp_osmr.xlsx", sheet = "one_adj_exp_憂鬱_out_IBS_top40")
format_snp_osmr_15 <- read_excel("format-snp_osmr.xlsx", sheet = "one_adj_exp_憂鬱_out_IBS_top15")

tp <- mtag.dep.two |>
  filter(str_detect(SNP, "^rs")) |>
  select(SNP, A1, A2, N, mtag_beta, mtag_se, mtag_pval, FRQ) 
write.table(tp, "mtag_gwas_snp/twosample/mtag/two-mtag-dep.txt", 
            sep = " ", quote = FALSE, row.names = FALSE)

tp <- gwas.ibs.two |>
  filter(str_detect(SNP, "^rs")) |>
  select(SNP, A1, NMISS, OR, SE, P) 
write.table(tp, "mtag_gwas_snp/twosample/gwas/two-gwas-ibs.txt", 
            sep = " ", quote = FALSE, row.names = FALSE)

# top
top_200 <- data.table(mtag.ibs)
top_200 <- top_200[order(mtag_pval)][1:200]
write_csv(top_200, 'mtag_gwas_snp/one-mtag-ibs.csv')

# fill in A2
dt <- data.table(gwas.dep)
get_A2 <- function(A1) {
  return(switch(A1, "A" = "T", "T" = "A", "C" = "G", "G" = "C"))
}
dt$A2 <- mapply(get_A2, dt$A1)

gwas.dep <- data.frame(dt)
gwas.dep <- apply(df, 2, as.character)
tp <- head(gwas.dep, 50)
write.table(tp, "mtag_gwas_snp/onesample/gwas/sample.txt", 
            sep = " ", quote = FALSE, row.names = FALSE)

# extract hg19 pos
library(rtracklayer)
library(GenomicRanges)
download.file("https://hgdownload.soe.ucsc.edu/gbdb/hg38/liftOver/hg38ToHg19.over.chain.gz", "hg38ToHg19.over.chain.gz")
system("gzip -d hg38ToHg19.over.chain.gz")

gwas_data <- gwas.dep
gr <- GRanges(seqnames = paste0("chr", gwas_data$CHR),
              ranges = IRanges(start = gwas_data$BP, end = gwas_data$BP),
              strand = "*",
              SNP = gwas_data$SNP)
chain <- import.chain("hg38ToHg19.over.chain")
gr_hg19 <- liftOver(gr, chain)
gr_hg19 <- unlist(gr_hg19)
gwas_data_hg19 <- data.frame(
  CHR = as.character(seqnames(gr_hg19)),
  BP_hg19 = start(gr_hg19),
  SNP = gr_hg19$SNP
)
gwas_data_hg19$CHR <- sub("^chr", "", gwas_data_hg19$CHR)

# ANNOVAR res
tp <- read.table('fuma/dep-one-m_524702/FUMA_job524702/snps.txt', header = TRUE)
dep.1m <- tp |> filter(!is.na(gwasP))
tp <- dep.1m |>
  group_by(nearestGene) |>
  filter(gwasP == min(gwasP, na.rm = TRUE)) |>
  ungroup()
write.table(tp, 'fuma/dep-1m.txt',
            sep = "\t", quote = FALSE, row.names = FALSE)

# pos + eQTL + ci
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