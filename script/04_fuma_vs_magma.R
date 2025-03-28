# see import of gwas.ded, gwas.dep, gwas.ibs in fuma format

library(tidyverse)

df <- inner_join(fuma_genes.ibs, magma_ibs, by = join_by('ensg' == 'GENE'))
write.table(df, 'fuma/fuma-magma_ibs.txt',
            sep = "\t", row.names = FALSE, quote = FALSE)