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

# ded vs. dep chord -------------------------------------------------------

# getting gene sets
ov <- read_excel('fuma/GO.xlsx', sheet = "GO_ded-dep") |> filter(Chord != 0)
gene_sets_of_interest <- ov$FULL_NAME

genes_in_sets <- c() # getting genes ensg
magma.genes_filtered <- magma_ded
for(set in gene_sets_of_interest){
  genes_list <- A_named_long[[set]]
  genes_filtered <- genes_list[genes_list %in% magma.genes_filtered$GENE]
  if(length(genes_filtered) == 0){
    gene_sets_of_interest <- gene_sets_of_interest[gene_sets_of_interest != set]
  }
  genes_in_sets <- c(genes_in_sets, genes_filtered)
}
genes_in_sets <- unique(genes_in_sets)
length(genes_in_sets)
length(gene_sets_of_interest) # 6
ded_genes <- genes_in_sets #75 .05 -> 21 .01 32
dep_genes <- genes_in_sets #83 .05 -> 24 .01 13
# no overlap
ov_genes <- unique(c(dep_genes, ded_genes)) # 154 .05 -> 45 .01 44

gnames <- c() # getting gene names
df <- filter(magma_dep, GENE %in% ov_genes) |> filter(!(GENE %in% ded_genes))
gnames_dep <- df$SYMBOL
df <- filter(magma_ded, GENE %in% ov_genes)
gnames_ded <- df$SYMBOL
gnames <- unique(c(gnames_dep, gnames_ded))
length(gnames)

# prepare matrix for graph
sub_matrix <- matrix(NA, nrow = length(ov_genes), ncol = length(gene_sets_of_interest),
                     dimnames = list(gnames, gene_sets_of_interest)) # r 45 * c 6

for (gene in rownames(sub_matrix)) {
  
  g <- -1
  p_dep <- magma_dep$P[magma_dep$SYMBOL == gene]
  if (length(p_dep) != 0) ensg <- magma_dep$GENE[magma_dep$SYMBOL == gene]
  
  p_ded <- magma_ded$P[magma_ded$SYMBOL == gene]
  if (length(p_ded) != 0) ensg <- magma_ded$GENE[magma_ded$SYMBOL == gene]
  
  if (length(p_ded) == 0) g <- -log(p_dep)
  if (length(p_dep) == 0) g <- -log(p_ded)
  if (g == -1) g <- (-log(p_dep) * -log(p_ded))
  
  for (gene_set in colnames(sub_matrix)) {
    
    p_dep_gs <- gsea.dep$P[gsea.dep$FULL_NAME == gene_set]
    p_ded_gs <- gsea.ded$P[gsea.ded$FULL_NAME == gene_set]
    gs <- -log(p_dep_gs) * -log(p_ded_gs)
    
    if(ensg %in% A_named_long[[gene_set]]) sub_matrix[gene, gene_set] <- g * gs
    else sub_matrix[gene, gene_set] <- 0
  }
}



# chord diagram
library(circlize)

# ibs-ded_go_5e-3_gene_1e-2
GO_label <- read.delim("script/GO_label")
mat <- sub_matrix
colnames(mat) <- GO_label$Names
gene_types <- list(Type1 = gnames_dep, Type2 = gnames_ded)

gene_set_color <- "grey"   
type1_color <- "#94C0FD" #DEP            
type2_color <- "#FE908F" #DED           

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

# ibs vs. ded chord -------------------------------------------------------

# getting gene sets
ov <- read_excel('fuma/GO.xlsx', sheet = "GO_ibs-ded") |> filter(Chord != 0)
gene_sets_of_interest <- ov$FULL_NAME 
gene_sets_of_interest <- gene_sets_of_interest[-2]

genes_in_sets <- c() # getting genes ensg
magma.genes_filtered <- magma_ded
for(set in gene_sets_of_interest){
  genes_list <- A_named_long[[set]]
  genes_filtered <- genes_list[genes_list %in% magma.genes_filtered$GENE]
  if(length(genes_filtered) == 0){
    gene_sets_of_interest <- gene_sets_of_interest[gene_sets_of_interest != set]
  }
  else {
    print(set)
    print(length(genes_filtered))
    if(length(genes_filtered) >= 5){
    genes_filtered <- magma.genes_filtered |>
      filter(GENE %in% genes_filtered) |>
      arrange(P) |>
      slice_head(n = 3) |>
      pull(GENE)
  }}
  genes_in_sets <- c(genes_in_sets, genes_filtered)
  # print(length(genes_in_sets))
}
genes_in_sets <- unique(genes_in_sets)
length(genes_in_sets)
length(gene_sets_of_interest) # 10
ded_genes <- genes_in_sets #78 .05 -> 21 .01 27
ibs_genes <- genes_in_sets #83 .05 -> 24 .01 23
# no overlap
ov_genes <- unique(c(ibs_genes, ded_genes)) # 154 .05 -> 45 .01 49
length(ov_genes)

gnames <- c() # getting gene names
df <- filter(magma_ibs, GENE %in% ibs_genes) |> filter(!(GENE %in% ded_genes))
gnames_ibs <- df$SYMBOL
df <- filter(magma_ded, GENE %in% ded_genes) 
gnames_ded <- df$SYMBOL
gnames <- unique(c(gnames_ibs, gnames_ded))
length(gnames)

# prepare matrix for graph
sub_matrix <- matrix(NA, nrow = length(ov_genes), ncol = length(gene_sets_of_interest),
                     dimnames = list(gnames, gene_sets_of_interest)) # r 45 * c 6

for (gene in rownames(sub_matrix)) {
  
  g <- -1
  p_ibs <- magma_ibs$P[magma_ibs$SYMBOL == gene]
  if (length(p_ibs) != 0) ensg <- magma_ibs$GENE[magma_ibs$SYMBOL == gene]

  p_ded <- magma_ded$P[magma_ded$SYMBOL == gene]
  if (length(p_ded) != 0) ensg <- magma_ded$GENE[magma_ded$SYMBOL == gene]

  if (length(p_ded) == 0) g <- -log(p_ibs)
  if (length(p_ibs) == 0) g <- -log(p_ded)
  if (g == -1) g <- (-log(p_ibs) * -log(p_ded))

  for (gene_set in colnames(sub_matrix)) {
    
    p_ibs_gs <- gsea.ibs$P[gsea.ibs$FULL_NAME == gene_set]
    p_ded_gs <- gsea.ded$P[gsea.ded$FULL_NAME == gene_set]
    gs <- -log(p_ibs_gs) * -log(p_ded_gs)

    if(ensg %in% A_named_long[[gene_set]]) sub_matrix[gene, gene_set] <- g * gs
    else sub_matrix[gene, gene_set] <- 0
  }
}



# chord diagram
library(circlize)

# ibs-ded_go_5e-3_gene_1e-2
GO_label <- c("Acute promyelocytic leukemia down", 
              "Production of molecular of immune response down", 
              # "Production of molecular mediator of immune response", 
              "Negative production of cytokine production involved in immune response", 
              "Regulation of neutrophil cheotaxis", 
              "Neural tube development", 
              "Regulation of production of molecular of immune response", 
              "Dendritic cell maturation up",        
              "MDM4 targets neuroepithelium down", 
              "SEMA4D-induced cell migration and growth cone collapse", 
              "MYC targets down")
mat <- sub_matrix
colnames(mat) <- GO_label
gene_types <- list(Type1 = gnames_ibs, Type2 = gnames_ded)

gene_set_color <- "grey"   
type1_color <- "#FEF197"            
type2_color <- "#FE908F"           

grid_colors <- c(rep(gene_set_color, ncol(mat)),  
                 rep(type1_color, length(gene_types$Type1)),  
                 rep(type2_color, length(gene_types$Type2)))

sector_order <- c(colnames(mat), gene_types$Type1, gene_types$Type2)

circos.clear()
circos.par(start.degree = 90)
chordDiagram(mat, annotationTrack = "grid", 
             preAllocateTracks = 1, grid.col = grid_colors, 
             order = sector_order)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = .375)
}, bg.border = NA)


# ibs vs. dep chord -------------------------------------------------------
ov <- read_excel('fuma/GO.xlsx', sheet = "GO_ibs-dep") |> filter(Chord != 0)
gene_sets_of_interest <- ov$FULL_NAME

genes_in_sets <- c() # getting genes ensg
magma.genes_filtered <- magma_ibs
for(set in gene_sets_of_interest){
  genes_list <- A_named_long[[set]]
  genes_filtered <- genes_list[genes_list %in% magma.genes_filtered$GENE]
  if(length(genes_filtered) == 0){
    gene_sets_of_interest <- gene_sets_of_interest[gene_sets_of_interest != set]
  }
  genes_in_sets <- c(genes_in_sets, genes_filtered)
}
genes_in_sets <- unique(genes_in_sets)
length(genes_in_sets)
length(gene_sets_of_interest) # 3
ibs_genes <- genes_in_sets #8
dep_genes <- genes_in_sets #11
ov_genes <- unique(c(ibs_genes, dep_genes)) #17

gnames <- c() # getting gene names
df <- filter(magma_ibs, GENE %in% ov_genes) |> filter(!(GENE %in% dep_genes))
gnames_ibs <- df$SYMBOL
df <- filter(magma_dep, GENE %in% ov_genes)
gnames_dep <- df$SYMBOL
gnames <- unique(c(gnames_ibs, gnames_dep))
length(gnames) # 17

# prepare matrix for graph
sub_matrix <- matrix(NA, nrow = length(ov_genes), ncol = length(gene_sets_of_interest),
                     dimnames = list(gnames, gene_sets_of_interest)) # r 26 * c 3

for (gene in rownames(sub_matrix)) {
  
  g <- -1
  p_ibs <- magma_ibs$P[magma_ibs$SYMBOL == gene]
  if (length(p_ibs) != 0) ensg <- magma_ibs$GENE[magma_ibs$SYMBOL == gene]
  
  p_dep <- magma_dep$P[magma_dep$SYMBOL == gene]
  if (length(p_dep) != 0) ensg <- magma_dep$GENE[magma_dep$SYMBOL == gene]
  
  if (length(p_dep) == 0) g <- -log(p_ibs)
  if (length(p_ibs) == 0) g <- -log(p_dep)
  if (g == -1) g <- (-log(p_ibs) * -log(p_dep))
  
  for (gene_set in colnames(sub_matrix)) {
    
    p_ibs_gs <- gsea.ibs$P[gsea.ibs$FULL_NAME == gene_set]
    p_dep_gs <- gsea.dep$P[gsea.dep$FULL_NAME == gene_set]
    gs <- -log(p_ibs_gs) * -log(p_dep_gs)
    
    if(ensg %in% A_named_long[[gene_set]]) sub_matrix[gene, gene_set] <- g * gs
    else sub_matrix[gene, gene_set] <- 0
  }
}

test <- "GOBP_NEGATIVE_REGULATION_OF_PROTEIN_SERINE_THREONINE_KINASE_ACTIVITY"
# debug: ENSG00000118997, DNAH7, ibs
# [1] "PATIL_LIVER_CANCER"                                                  
# [2] "REACTOME_DEUBIQUITINATION"                                           
# [3] "GOBP_NEGATIVE_REGULATION_OF_PROTEIN_SERINE_THREONINE_KINASE_ACTIVITY"

# chord diagram
library(circlize)

# ibs-dep_go_5e-3_gene_1e-2
GO_label <- c("Scar complex", 
              "Positive regulation of APR2/3 complex mediated actin nucleation", 
              "Retinoblastoma")
mat <- sub_matrix
colnames(mat) <- GO_label
gene_types <- list(Type1 = gnames_ibs, Type2 = gnames_dep)

gene_set_color <- "grey"   
type1_color <- "#FEF197"  # IBS          
type2_color <- "#94C0FD" # DEP          

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
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = .7)
}, bg.border = NA)

circos.clear()