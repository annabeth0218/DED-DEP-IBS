library(ggplot2)
library(dplyr)
library(readr)
library(ggtext)

# select map
df <- read_table('figs/eqtl_heatmap/eqtlmaps.txt') |>
  filter(relevance == 2 | relevance ==3) |>
  filter(is.na(GTEx) | GTEx == 1)
write.table(df, 'figs/eqtl_heatmap/eqtlmaps_used.txt',
            quote = F, row.names = F)

# get gene list, map list
df <- read_table('figs/eqtl_heatmap/genes.txt')
genes <- df$gene

df <- read_table('figs/eqtl_heatmap/eqtlmaps_used.txt')
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
eqtl_ded <- process_eqtl_data("fuma/fuma_ded/eqtl.txt", "#C14F4F")  # Dry Eye (DED)
eqtl_dep <- process_eqtl_data("fuma/fuma_dep/eqtl.txt", "#4A5680") # Depression (DEP)
eqtl_ibs <- process_eqtl_data("fuma/fuma_ibs/eqtl.txt", "#D19A4F") # Irritable Bowel Syndrome (IBS)

# Combine data
combined_data <- bind_rows(eqtl_ded, eqtl_dep, eqtl_ibs)

# Convert categorical variables to factors for plotting
combined_data$Category <- factor(combined_data$Category, levels = c("#C14F4F", "#4A5680", "#D19A4F"))

# Ensure `!is.na(RiskIncAllele)` dots are plotted last (on top)
combined_data <- combined_data |> arrange(not_risk_allele)
combined_data[combined_data$tissue == 'CRBL', 'signed_stats'] <- 0.1
# combined_data <- combined_data |> filter(!(symbol == 'EYS' & eqtlmap == 'BRAINEAC#CRBL'))
combined_data$symbol_colored <- paste0("<span style='color:", combined_data$Category, "'>", combined_data$symbol, "</span>")
combined_data$symbol_colored <- factor(combined_data$symbol_colored, levels = unique(combined_data$symbol_colored[order(match(combined_data$symbol, genes))]))

p <- ggplot(combined_data, aes(x = factor(eqtlmap, levels = eqtlmaps), 
                               y = symbol_colored, 
                               size = -log10(p),
                               fill = signed_stats,
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
    axis.text.y = element_markdown(vjust = 0.5, size = 16, face = "bold"),  # Y-axis with colo#B57979 labels
    axis.ticks = element_blank()
  ) +
  labs(fill = "Signed Stats", size = "-log10(P-value)", x = NULL, y = NULL)

# Save plot
ggsave("figs/eqtl_heatmap/heatmap_eqtl.png", plot = p, width = 24, height = 14, dpi = 300)
