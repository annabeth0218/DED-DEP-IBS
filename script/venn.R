# tutorial https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/
library(ggvenn)
library(data.table)

venn2 <- list(
  FUMA = fuma_genes.ibs$ensg,
  MAGMA = magma_ibs$GENE
)

ggvenn(
  venn2, 
  fill_color = c("#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 5
)

venn3 <- list(
  DED = fuma_genes.ded$ensg,
  DEP = fuma_genes.dep$ensg,
  IBS = fuma_genes.ibs$ensg
)

# color: "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"
ggvenn(
  venn3, 
  fill_color = c("#0073C2FF", "#CD534CFF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 5
)



