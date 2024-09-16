# tutorial https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/
library(ggvenn)
library(data.table)

top_200 <- data.table(mtag.ibs)
t.mtag.ibs <- top_200[order(mtag_pval)][1:200]

mg <- list(
  MTAG <- t.mtag.ded$SNP,
  GWAS <- t.mtag.ded$SNP,
  MTAG2 <- t.mtag.ded$SNP
)

# color: "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"
ggvenn(
  mg, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 5
)


