library(tidyverse)

df <- read_table('fuma/fuma_ded/FUMA_celltype585207/magma_celltype_step1.txt')
celltype <- df

celltype <- celltype |>
  filter(P < 0.05)