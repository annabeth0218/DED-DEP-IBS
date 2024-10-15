library(tidyverse)

gsea.ibs <- read_excel("fuma/magma.xlsx", skip = 4, sheet = "ibs-1g.gsa.out") #16987
gsea.ded <- read_excel("fuma/magma.xlsx", skip = 4, sheet = "ded-1g.gsa.out") #16986
gsea.dep <- read_excel("fuma/magma.xlsx", skip = 4, sheet = "dep-1g.gsa.out") #16986

gsea.ded <- gsea.ded |> filter(P < 0.05) #874
gsea.dep <- gsea.dep |> filter(P < 0.05) #976
gsea.ibs <- gsea.ibs |> filter(P < 0.05) #863

ov.ded.dep <- semi_join(gsea.ded, gsea.dep, by = "FULL_NAME")
ov.ibs.dep <- semi_join(gsea.ibs, gsea.dep, by = "FULL_NAME")
ov.ibs.ded <- semi_join(gsea.ibs, gsea.ded, by = "FULL_NAME")
write.table(ov, "fuma/GO_ibs-dep.txt", quote = FALSE, row.names = FALSE)

msigdb <- read_excel("fuma/MSigDB_20231Hs_MAGMA.xlsx", col_names = TRUE)
colnames(msigdb) <- c("GeneSet", "UnwantedInfo", paste0("Gene", 1:(ncol(msigdb)-2)))
A <- msigdb |> select(-UnwantedInfo)
A_list <- apply(A, 1, function(row) as.list(row[!is.na(row)]))

ov.ded.dep$GENES <- sapply(ov.ded.dep$FULL_NAME, function(gene_set) {
  set_in_A <- A_list[sapply(A_list, function(x) x[[1]] == gene_set)]
  if (length(set_in_A) == 1) {
    genes_in_set <- unlist(set_in_A)[-1] 
    overlapping_genes <- intersect(genes.ded.1g$ensg, genes_in_set)
    return(paste(overlapping_genes, collapse = ", "))
  } else {
    return(NA)  
  }
})
  

# overlap
overlap <- anti_join(gwas.ded, gwas.dep, by = "SNP") |>
  filter(P < 0.05)

snps.ded.1g <- read_excel("fuma/snp2gene.xlsx", sheet = "ded-1g.snps")
overlap_updated <- left_join(overlap, snps.ded.1g, by = c("SNP" = "rsID")) 
overlap_updated <- overlap_updated |> filter(!is.na(uniqID))


