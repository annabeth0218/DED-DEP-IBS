# install.packages("BiocManager")
BiocManager::install("snpStats")
install.packages("remotes")
remotes::install_github("josefin-werme/LAVA")
library(LAVA)

# A2
dt <- data.table(gwas.ibs)
get_A2 <- function(A1) {
  return(switch(A1, "A" = "T", "T" = "A", "C" = "G", "G" = "C"))
}
dt$A2 <- mapply(get_A2, dt$A1)

df <- data.frame(dt)
df <- apply(df, 2, as.character)

write.table(tp, "mtag_gwas_snp/onesample/gwas/sample.txt", 
            sep = " ", quote = FALSE, row.names = FALSE)

out <- df |>
  select(SNP, A1, A2, NMISS, STAT)
write.table()