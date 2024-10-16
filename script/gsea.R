library(tidyverse)

# import GO results of FUMA
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
ov.gs <- unique(c(ov.ded.dep$FULL_NAME, ov.ibs.ded$FULL_NAME, ov.ibs.dep$FULL_NAME)) #159
gs <- unique(c(gsea.ded$FULL_NAME, gsea.dep$FULL_NAME, gsea.ibs$FULL_NAME)) #2551

# complete gene lists of all gene sets
msigdb <- read_excel("fuma/MSigDB_20231Hs_MAGMA.xlsx", col_names = TRUE) #17022
colnames(msigdb) <- c("GeneSet", "UnwantedInfo", paste0("Gene", 1:(ncol(msigdb)-2)))
A <- msigdb |>
  select(-UnwantedInfo) |>
  filter(GeneSet %in% gs)
A_long <- lapply(1:nrow(A), function(i) {
  genes <- as.list(A[i, -1][!is.na(A[i, -1])])
  setNames(list(genes), A$GeneSet[i])
})
A_named_long <- do.call(c, A_long)

# ov <- ov.ibs.ded |> #53
#   filter(P < 0.05) |>
#   filter(NGENES < 40) #10
gene_sets_of_interest <- ov$FULL_NAME

length(A_named_long[['GOMF_STEROID_DEHYDROGENASE_ACTIVITY_ACTING_ON_THE_CH_OH_GROUP_OF_DONORS_NAD_OR_NADP_AS_ACCEPTOR']]) 
# 31 <- 24

genes_in_sets <- c()
for(set in gene_sets_of_interest){
  genes_list <- A_named_long[[set]]
  genes_filtered <- genes_list[genes_list %in% magma_ibs$GENE]
  genes_in_sets <- c(genes_in_sets, genes_filtered)
}
genes_in_sets <- unique(genes_in_sets)

test <- magma_ibs |> # 3
  filter(GENE %in% A_named_long[['GOMF_STEROID_DEHYDROGENASE_ACTIVITY_ACTING_ON_THE_CH_OH_GROUP_OF_DONORS_NAD_OR_NADP_AS_ACCEPTOR']])

testgs <- 'GOMF_STEROID_DEHYDROGENASE_ACTIVITY_ACTING_ON_THE_CH_OH_GROUP_OF_DONORS_NAD_OR_NADP_AS_ACCEPTOR'

genes_in_sets_p <- genes_in_sets |> filter()

sub_matrix <- matrix(NA, nrow = length(genes_in_sets), ncol = length(gene_sets_of_interest),
                     dimnames = list(genes_in_sets, gene_sets_of_interest))

sub_matrix <- gene_matrix[rownames(gene_matrix) %in% gene_sets_of_interest, ]
genes_in_sets <- lapply(rownames(sub_matrix), function(set) {
  present_genes <- colnames(sub_matrix)[sub_matrix[set, ] == 1]
  return(present_genes)
})
names(genes_in_sets) <- rownames(sub_matrix)



# chord diagram
library(circlize)

# saneky diagram

# overlap
overlap <- anti_join(gwas.ded, gwas.dep, by = "SNP") |>
  filter(P < 0.05)

snps.ded.1g <- read_excel("fuma/snp2gene.xlsx", sheet = "ded-1g.snps")
overlap_updated <- left_join(overlap, snps.ded.1g, by = c("SNP" = "rsID")) 
overlap_updated <- overlap_updated |> filter(!is.na(uniqID))

ggplot(ov., aes(x = FULL_NAME, y = P)) +
  geom_point() +
  theme(axis.text.x = element_blank())
