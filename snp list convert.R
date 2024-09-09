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

write.table(gwas.ded, "mtag_gwas_snp/onesample/gwas/gwas_ded.txt", 
            sep = " ", quote = FALSE, row.names = FALSE)