library(tidyverse)

# original hg38 build
# gwas.dep <- read.csv('table/gwas_snp/with_A2/one_sample_gwas/depression_log.csv')
# gwas.ded <- read.csv('table/gwas_snp/with_A2/one_sample_gwas/dryeye_case60_log.csv')
# gwas.ibs <- read.csv('table/gwas_snp/with_A2/one_sample_gwas/ibs_case50_log.csv')

gwas.dep <- read.table('table/gwas_snp/with_A2/one_sample_gwas/hg19_gwas_dep.txt')
gwas.ded <- read.table('table/gwas_snp/with_A2/one_sample_gwas/hg19_gwas_ded.txt')
gwas.ibs <- read.table('table/gwas_snp/with_A2/one_sample_gwas/hg19_gwas_ibs.txt')
