library(tidyverse)

# original hg38 build
# gwas.dep <- read.csv('table/gwas_snp/with_A2/one_sample_gwas/depression_log.csv')
# gwas.ded <- read.csv('table/gwas_snp/with_A2/one_sample_gwas/dryeye_case60_log.csv')
# gwas.ibs <- read.csv('table/gwas_snp/with_A2/one_sample_gwas/ibs_case50_log.csv')

gwas.dep <- read.table('table/gwas_snp/with_A2/one_sample_gwas/hg19_gwas_dep.txt', header = T)
gwas.ded <- read.table('table/gwas_snp/with_A2/one_sample_gwas/hg19_gwas_ded.txt', header = T)
gwas.ibs <- read.table('table/gwas_snp/with_A2/one_sample_gwas/hg19_gwas_ibs.txt', header = T)

magma_ibs <- read.delim("fuma/fuma_ibs/magma.genes.out") |> filter(P < 0.05)
magma_ded <- read.delim("fuma/fuma_ded/magma.genes.out") |> filter(P < 0.05)
magma_dep <- read.delim("fuma/fuma_dep/magma.genes.out") |> filter(P < 0.05)

snps.ded <- read.table("fuma/fuma_ded/snps.txt", header = T)
snps.dep <- read.table("fuma/fuma_dep/snps.txt", header = T)
snps.ibs <- read.table("fuma/fuma_ibs/snps.txt", header = T)

riskloci.ded <- read.table("fuma/fuma_ded/IndSigSNPs.txt", header = T)
riskloci.dep <- read.table("fuma/fuma_dep/IndSigSNPs.txt", header = T)
riskloci.ibs <- read.table("fuma/fuma_ibs/IndSigSNPs.txt", header = T)

df <- snps.ded[snps.ded$rsID %in% riskloci.ded$rsID, ]
write.table(df, file = "table/Riskloci_ded.txt", sep = "\t", row.names = FALSE, quote = FALSE)
df <- snps.dep[snps.dep$rsID %in% riskloci.dep$rsID, ]
write.table(df, file = "table/Riskloci_dep.txt", sep = "\t", row.names = FALSE, quote = FALSE)
df <- snps.ibs[snps.ibs$rsID %in% riskloci.ibs$rsID, ]
write.table(df, file = "table/Riskloci_ibs.txt", sep = "\t", row.names = FALSE, quote = FALSE)
