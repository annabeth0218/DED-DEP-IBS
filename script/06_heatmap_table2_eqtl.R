library(ggplot2)
library(dplyr)
library(readr)
library(ggtext)

df <- read_table('figs/fig2_eqtl_heatmap/genes.txt')
genes <- df$gene
# df <- read_table('figs/fig2_eqtl_heatmap/eqtlmaps.txt') |>
#   filter(relevance == 2 | relevance ==3) |>
#   filter(is.na(GTEx) | GTEx == 1)
# write.table(df, 'figs/fig2_eqtl_heatmap/eqtlmaps_used.txt',
#             quote = F, row.names = F)
df <- read_table('figs/fig2_eqtl_heatmap/eqtlmaps_v4.txt')
eqtlmaps <- df$eqtlmap

# Function to process each file
process_eqtl_data <- function(file_path, category_color) {
  df <- read_tsv(file_path) |>
    mutate(eqtlmap = paste(db, tissue, sep = '#')) |>
    filter(eqtlmap %in% eqtlmaps) |>
    filter(symbol %in% genes) |> 
    mutate(not_risk_allele = !is.na(RiskIncAllele))  # Create flag for visibility
  
  df_filtered <- df |>
    group_by(gene, eqtlmap) |>
    slice_min(order_by = p, n = 1) |>
    ungroup()
  
  df_filtered <- df_filtered |> mutate(Category = category_color)
  return(df_filtered)
}

# Load data
eqtl_ded <- process_eqtl_data("fuma/fuma_ded/eqtl.txt", "#BF0000")  # Dry Eye (DED)
eqtl_dep <- process_eqtl_data("fuma/fuma_dep/eqtl.txt", "#576FA0") # Depression (DEP)
eqtl_ibs <- process_eqtl_data("fuma/fuma_ibs/eqtl.txt", "#BA8E23") # Irritable Bowel Syndrome (IBS)

# Combine data
combined_data <- bind_rows(eqtl_ded, eqtl_dep, eqtl_ibs)

# Convert categorical variables to factors for plotting
combined_data$Category <- factor(combined_data$Category, levels = c("#BF0000", "#576FA0", "#BA8E23"))

# Ensure `!is.na(RiskIncAllele)` dots are plotted last (on top)
combined_data <- combined_data |> arrange(not_risk_allele)
combined_data[combined_data$tissue == 'CRBL', 'signed_stats'] <- -0.337
# combined_data <- combined_data |> filter(!(symbol == 'EYS' & eqtlmap == 'BRAINEAC#CRBL'))
combined_data$symbol_colored <- paste0("<span style='color:", combined_data$Category, "'>", combined_data$symbol, "</span>")
combined_data$symbol_colored <- factor(combined_data$symbol_colored, levels = unique(combined_data$symbol_colored[order(match(combined_data$symbol, genes))]))
names(combined_data)[names(combined_data) == "signed_stats"] <- "Beta"

p <- ggplot(combined_data, aes(x = factor(eqtlmap, levels = eqtlmaps), 
                               y = symbol_colored, 
                               size = -log10(p),
                               fill = Beta,
                               color = ifelse(not_risk_allele, "black", "grey"))) + 
  geom_point(shape = 21, 
             stroke = ifelse(combined_data$not_risk_allele, 1.2, 0.5)) + 
  scale_size(range = c(4, 8)) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  
  scale_color_identity() +  
  scale_x_discrete(position = "top") +  # Move x-axis labels to the top
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = "grey80"),  
    plot.background = element_rect(fill = "white", color = NA),  
    axis.line = element_blank(),  
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 1),  
    legend.key = element_blank(),  
    legend.background = element_blank(),  
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, 
                                   size = 16, face = "bold", color = "#363737"), # X-axis at top
    axis.text.y = element_markdown(vjust = 0.5, size = 16, face = "bold"),  # Y-axis with colored labels
    axis.ticks = element_blank()
  ) +
  labs(fill = "Beta", size = "-log10(P-value)", x = NULL, y = NULL)

# Save plot
# all: 24*14
# part: 14*14
ggsave("figs/fig2_eqtl_heatmap/fig2_.png", plot = p, width = 14, height = 14, dpi = 300)


# EYS issue
# tutorial https://adairama.wordpress.com/2013/10/25/278/
load("figs/fig2_eqtl_heatmap/eqtlmaps/2959403.rda")
length(expr)
table(a <- colMeans(is.na(expr[["CRBL"]])))
length(valids <- names( which(a==0) ))
library(MatrixEQTL)
eys.expr <- SlicedData$new()
eys.expr$CreateFromMatrix( as.matrix( expr[["CRBL"]][ , valids ] ) )

eys.markers <- SlicedData$new()
eys.markers$CreateFromMatrix( as.matrix( markers[ , valids ] ) )
identical( colnames(eys.expr), colnames(eys.markers) ) #TRUE

tmpf <- tempfile() #create a place to write output file

res <- Matrix_eQTL_main( eys.markers, eys.expr, cvrt=SlicedData$new(),
                            output_file_name=tmpf, pvOutputThreshold=1e-5,
                            useModel=modelLINEAR, errorCovariance=numeric(0) )
unlink(tmpf) #del the tmp dir
out <- res$all$eqtls
out <- out[ , c("snps", "gene", "beta", "statistic", "pvalue")]
dim(out)
#from eqtl_ibs, search 6:66169835:C:T, get beta=-0.3370683

