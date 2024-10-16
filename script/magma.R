library(tidyverse)

# for ibs-1g r^2 6e-1
magma_ibs <- read.delim("~/R_Projects/ded-dep-ibs/fuma/r2_6e-1/ibs-1g_524698/FUMA_job524698/magma.genes.out") |>
  filter(P < 0.05)
nrow(magma_ibs) #15763 -> 930
# for ded-1g r^2 6e-1
magma_ded <- read.delim("~/R_Projects/ded-dep-ibs/fuma/r2_6e-1/ded-1g_524703/FUMA_job524703/magma.genes.out") |>
  filter(P < 0.05)
nrow(magma_ded) #15763 -> 1120
# for dep-1g r^2 6e-1
magma_dep <- read.delim("~/R_Projects/ded-dep-ibs/fuma/r2_6e-1/dep-1g_524690/FUMA_job524690/magma.genes.out") |>
  filter(P < 0.05)
nrow(magma_dep) #15759 -> 485

magma.genes_filtered <- magma.genes |>
  filter(P < 0.05) |>
  filter(NSNPS >= 3)
nrow(magma.genes_filtered) #1053, 981, 940

ov <- semi_join(magma.genes_filtered, genes.ded.1g, by = c("GENE" = "ensg"))

# MUC16 ENSG00000181143
magma_ded |>
  filter(GENE == "ENSG00000181143") #found
magma_dep |>
  filter(GENE == "ENSG00000181143") #no, P 0.055965
magma_ibs |>
  filter(GENE == "ENSG00000181143") #no, P 0.84454

# TENM2 ENSG00000145934
magma_ded |>
  filter(GENE == "ENSG00000145934") #found
magma_dep |>
  filter(GENE == "ENSG00000145934") #no, 0.83197
magma_ibs |>
  filter(GENE == "ENSG00000145934") #no, 0.82849

# venn
hvh <- list(
  # DED = magma_ded$GENE,
  MAGMA = magma_ibs$GENE,
  FUMA = genes.ibs.1g$ensg
)

# color: "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"
ggvenn(
  hvh, 
  fill_color = c("#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 5
)