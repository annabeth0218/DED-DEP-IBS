library(data.table)
library(readxl)
library(BSgenome.Hsapiens.UCSC.hg38) # using GRCh38/hg38

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

# top
top_200 <- data.table(mtag.ibs)
top_200 <- top_200[order(mtag_pval)][1:200]
write_csv(top_200, 'mtag_gwas_snp/one-mtag-ibs.csv')

# prep for FUMA
dt <- data.table(gwas.dep)
get_A2 <- function(A1) {
  return(switch(A1, "A" = "T", "T" = "A", "C" = "G", "G" = "C"))
}
dt$A2 <- mapply(get_A2, dt$A1)

gwas.dep <- data.frame(dt)
gwas.dep <- apply(df, 2, as.character)
write.table(gwas.dep, "mtag_gwas_snp/onesample/gwas/gwas_dep.txt", 
            sep = " ", quote = FALSE, row.names = FALSE)


