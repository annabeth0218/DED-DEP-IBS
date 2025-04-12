# eda
# ov.ded.dep <- semi_join(gsea.ded, gsea.dep, by = "FULL_NAME")
# ov.ibs.dep <- semi_join(gsea.ibs, gsea.dep, by = "FULL_NAME")
# ov.ibs.ded <- semi_join(gsea.ibs, gsea.ded, by = "FULL_NAME")
# write.table(ov, "fuma/GO_ibs-dep.txt", quote = FALSE, row.names = FALSE)
# ov.gs <- unique(c(ov.ded.dep$FULL_NAME, ov.ibs.ded$FULL_NAME, ov.ibs.dep$FULL_NAME)) #159
# gs <- unique(c(gsea.ded$FULL_NAME, gsea.dep$FULL_NAME, gsea.ibs$FULL_NAME)) #2551

library(tidyverse)
library(readxl)
library(circlize)

# import GO results of FUMA
gsea.ibs <- read_excel("fuma/magma_GO.xlsx", skip = 4, sheet = "ibs-1g.gsa.out") |> #16987
  filter(P < 0.05) #874
gsea.ded <- read_excel("fuma/magma_GO.xlsx", skip = 4, sheet = "ded-1g.gsa.out") |> #16986
  filter(P < 0.05) #976
gsea.dep <- read_excel("fuma/magma_GO.xlsx", skip = 4, sheet = "dep-1g.gsa.out") |> #16986
  filter(P < 0.05) #863

color_dep <- "#576FA0"             
color_ded <- "#B57979"
color_ibs <- "#E3B87F"

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

chord_get_genes <- function(magma_data, gene_sets, gene_mapping, run) {
  genes_in_sets <- c()
  filtered_gene_sets <- gene_sets
  
  for(set in gene_sets) {
    genes_list <- gene_mapping[[set]]
    genes_filtered <- genes_list[genes_list %in% magma_data$GENE]
    
    if(length(genes_filtered) == 0) {
      filtered_gene_sets <- filtered_gene_sets[filtered_gene_sets != set]
    }
    else if(run){
      if(length(genes_filtered) >= 5){
        genes_filtered <- magma_data |>
          filter(GENE %in% genes_filtered) |>
          arrange(P) |>
          slice_head(n = 3) |>
          pull(GENE)
      }}
    genes_in_sets <- c(genes_in_sets, genes_filtered)
  }
  list(genes = unique(genes_in_sets), genesets = filtered_gene_sets)
}

chord <- function(gs, df1, df2, gsea1, gsea2, type1_color, type2_color, run) {
  
  circos.clear()
  gene_sets_of_interest <- gs
  
  # getting genes
  res <- chord_get_genes(df1, gene_sets_of_interest, A_named_long, run) 
  genes1 <- res$genes
  gene_sets_of_interest <- res$genesets
  res <- chord_get_genes(df2, gene_sets_of_interest, A_named_long, run) 
  genes2 <- res$genes
  gene_sets_of_interest <- res$genesets
  res <- chord_get_genes(df1, gene_sets_of_interest, A_named_long, run) 
  genes1 <- res$genes
  gene_sets_of_interest <- res$genesets
  ov_genes <- unique(c(genes1, genes2)) 
  
  # getting gene names
  gnames <- c() 
  df <- filter(df1, GENE %in% genes1) |> filter(!(GENE %in% genes2))
  gnames1 <- df$SYMBOL
  df <- filter(df2, GENE %in% genes2)
  gnames2 <- df$SYMBOL
  gnames <- unique(c(gnames1, gnames2))
  
  # cat(length(ov_genes), '\n')
  # cat(length(gnames), '\n')
  # cat(length(gnames1), '\n')
  # cat(length(gnames2), '\n')
  # cat(length(genes1), '\n')
  # cat(length(genes2), '\n')
  # cat(length(gene_sets_of_interest), '\n')
  

  # prepare matrix for graph
  mat <- matrix(NA, nrow = length(ov_genes), ncol = length(gene_sets_of_interest),
                dimnames = list(gnames, gene_sets_of_interest)) 
  
  for (gene in rownames(mat)) {
    
    g <- -1
    p1 <- df1$P[df1$SYMBOL == gene]
    if (length(p1) != 0) ensg <- df1$GENE[df1$SYMBOL == gene]
    
    p2 <- df2$P[df2$SYMBOL == gene]
    if (length(p2) != 0) ensg <- df2$GENE[df2$SYMBOL == gene]
    
    if (length(p1) == 0) g <- -log(p2)
    if (length(p2) == 0) g <- -log(p1)
    if (g == -1) g <- (-log(p1) * -log(p2))
    
    for (gene_set in colnames(mat)) {
      
      p1_gs <- gsea1$P[gsea1$FULL_NAME == gene_set]
      p2_gs <- gsea2$P[gsea2$FULL_NAME == gene_set]
      gs <- -log(p1_gs) * -log(p2_gs)
      
      if(ensg %in% A_named_long[[gene_set]]) mat[gene, gene_set] <- g * gs
      else mat[gene, gene_set] <- 0
    }
  }
  
  # chord diagram
  colnames(mat) <- GO_label$Names
  gene_types <- list(Type1 = gnames1, Type2 = gnames2)
  
  gene_set_color <- "grey"   
  
  grid_colors <- c(rep(gene_set_color, ncol(mat)),  
                   rep(type1_color, length(gene_types$Type1)),  
                   rep(type2_color, length(gene_types$Type2)))
  
  sector_order <- c(colnames(mat), gene_types$Type1, gene_types$Type2)
  
  circos.par(start.degree = 90)
  chordDiagram(mat, annotationTrack = "grid", 
               preAllocateTracks = 1, grid.col = grid_colors, 
               order = sector_order)
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = .4)
  }, bg.border = NA)
  
}

# ded 75 .05 -> 21 .01 32
# dep 83 .05 -> 24 .01 13
# total genes 154 .05 -> 45 .01 44

ov <- read_excel('fuma/magma_GO.xlsx', sheet = "GO_deddep") |> filter(Chord != 0)
GO_label <- read.delim("script/GO_label") 
chord(ov$FULL_NAME, magma_ded, magma_dep, gsea.ded, gsea.dep, color_ded, color_dep, F)

# ded 78 .05 -> 21 .01 27
# ibs 83 .05 -> 24 .01 23
# total genes 154 .05 -> 45 .01 49

ov <- read_excel('fuma/magma_GO.xlsx', sheet = "GO_ibsded") |> filter(Chord != 0)
gs <- ov$FULL_NAME
gs <- gs[c(-2,-3)] # remove 2nd one
GO_label <- read.delim("script/GO_label") 
chord(gs, magma_ibs, magma_ded, gsea.ibs, gsea.ded, color_ibs, color_ded, T)

# ibs 8
# dep 11
# total genes 17

ov <- read_excel('fuma/magma_GO.xlsx', sheet = "GO_ibsdep") |> filter(Chord != 0)
GO_label <- read.delim("script/GO_label") 
chord(ov$FULL_NAME, magma_ibs, magma_dep, gsea.ibs, gsea.dep, color_ibs, color_dep, F)
