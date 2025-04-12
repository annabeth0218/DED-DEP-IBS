# tutorial https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/
library(ggvenn)
library(data.table)

ibs <- read.table("/Volumes/256_Anna/3MR/ibs_case50_new_QC2.bim", header = FALSE)
ded <- read.table("/Volumes/256_Anna/3MR/dryeye_case60_new_QC2.bim", header = FALSE)
dep <- read.table("/Volumes/256_Anna/3MR/depression_new_QC2.bim", header = FALSE)


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
  DED = ded$V2,
  DEP = dep$V2,
  IBS = ibs$V2
)

# color: "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"
ggvenn(
  venn3, 
  fill_color = c("#0073C2FF", "#CD534CFF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 5, text_size = 2
)



