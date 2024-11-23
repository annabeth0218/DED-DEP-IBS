# tutorial https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/
library(ggvenn)
library(data.table)

top_200 <- data.table(mtag.ibs)
t.mtag.ibs <- top_200[order(mtag_pval)][1:200]

mg200 <- list(
  MTAG = f.dep.1m$ensg,
  GWAS = f.dep.1g$ensg
)

tvt <- list(
  IBS = ibs_genes,
  DED = ded_genes
)

ggvenn(
  tvt, 
  fill_color = c("#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 5
)

hvh <- list(
  DED = magma_ded$GENE,
  DEP = magma_dep$GENE,
  IBS = magma_ibs$GENE
)

# color: "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"
ggvenn(
  hvh, 
  fill_color = c("#0073C2FF", "#CD534CFF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 5
)


