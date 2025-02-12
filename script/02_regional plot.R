library(ggplot2)
library(dplyr)
library(ggnewscale)

# Load your GWAS data (modify the file path and structure as needed)
# Example CSV file structure should include: CHR, BP, SNP, P
gwas_data <- read.csv("mtag_gwas/onesample/gwas/depression_log_ADD.csv")
g2 <- read.csv("mtag_gwas/onesample/gwas/ibs_case50_log_ADD.csv")
g3 <- read.csv("mtag_gwas/onesample/gwas/dryeye_case60_log_ADD.csv")

# Define the region to plot (modify these values as needed)
target_chr <- 10         # Example chromosome
start_pos <- 19890544 #20390544 #20679473   # Start base pair position
end_pos <- 23867387 #23367387 #23656316     # End base pair position

# Filter data to include only the target region
region_data <- gwas_data %>%
  filter(CHR == target_chr, BP >= start_pos, BP <= end_pos)
r2 <- g2 %>%
  filter(CHR == target_chr, BP >= start_pos, BP <= end_pos)
r3 <- g3 %>%
  filter(CHR == target_chr, BP >= start_pos, BP <= end_pos)

# Calculate -log10(p-value) for better visualization
region_data <- region_data %>%
  mutate(logP = -log10(P))

r2 <- r2 %>%
  mutate(logP = -log10(P))

r3 <- r3 %>%
  mutate(logP = -log10(P))

# Merge data to identify common SNPs and apply the p-value condition
merged_data <- region_data %>%
  inner_join(r2, by = c("CHR", "BP", "SNP"), suffix = c("_gwas1", "_gwas2")) %>%
  mutate(
    keep_in_gwas1 = P_gwas1 <= P_gwas2 & P_gwas1 > 1e-3 & P_gwas2 > 1e-3,
    keep_in_gwas2 = P_gwas2 < P_gwas1 & P_gwas2 > 1e-3 & P_gwas1 > 1e-3
  )

# Remove overlapping SNPs based on the condition
region_data <- region_data %>%
  filter(!(SNP %in% merged_data$SNP & merged_data$keep_in_gwas2))

r2 <- r2 %>%
  filter(!(SNP %in% merged_data$SNP & merged_data$keep_in_gwas1))

# Calculate -log10(p-value) for better visualization
region_data <- region_data %>%
  mutate(logP = -log10(P))

r2 <- r2 %>%
  mutate(logP = -log10(P))

# Visualization
ggplot() +
  # First GWAS data with red gradient
  geom_point(
    data = region_data, 
    aes(x = BP, y = logP, color = logP), 
    size = 2
  ) +
  scale_color_gradient(
    low = "darkgray", high = "red",
    name = "DEP"
  ) +
  
  # Add a new color scale for the second GWAS data
  new_scale_color() +
  
  # Second GWAS data with blue gradient
  geom_point(
    data = r2, 
    aes(x = BP, y = logP, color = logP), 
    size = 2, alpha = 0.7
  ) +
  scale_color_gradient(
    low = "darkgray", high = "blue",
    name = "IBS"
  ) +
  
  new_scale_color() +
  
  # Third GWAS data with blue gradient
  geom_point(
    data = r3, 
    aes(x = BP, y = logP, color = logP), 
    size = 2, alpha = 0.7
  ) +
  scale_color_gradient(
    low = "darkgray", high = "yellow",
    name = "DED"
  ) +
  
  # Labels and theme settings
  theme_minimal() +
  labs(
    x = "Base Pair Position",
    y = "-log10 P"
  ) +
  theme(
    legend.position = "top",
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  ) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(
    breaks = seq(0, 6.5, by = 0.5),
    limits = c(0, 6.5),
    expand = c(0, 0))