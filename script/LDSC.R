
gwas.ded |>
  select(SNP, A1, A2, NMISS, P, STAT) |>
  write.table('/Users/annabethlu/ldsc/data/sumstats_ded.txt', 
              sep = '\t', row.names = F, quote = F)