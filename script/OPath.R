# 1108 abstract

library(ggplot2)
library(dplyr)
library(tidyr)

# Data
metrics <- data.frame(
  Metric = c("BLEU1", "BLEU4", "ROUGE", "CIDEr-R", "METEOR"),
  Mean = c(0.3785, 0.2151, 0.4210, 0.0953, 0.5435),
  Min = c(0.2258, 0.0000, 0.2885, 0.0953, 0.3204),
  Max = c(0.5192, 0.3523, 0.5594, 0.0953, 0.7414)
)

# Plot
metrics$Metric <- factor(metrics$Metric, levels = unique(metrics$Metric))
metric_colors <- c(
  "BLEU1" = "#2AB9B3",  # #436b94
  "BLEU4" = "#69CEC9", # #829cb8
  "ROUGE" = "gray", # #ee6b6e
  "CIDEr-R" = "#e0a5a6", # #e8b028
  "METEOR" = "#B10002" # #6a4a3a
)

ggplot(metrics, aes(x = Metric, y = Mean, fill = Metric)) + 
  geom_col(width = 1, color = "black") +  
  # geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, color = "gray") +
  geom_text(aes(label = sprintf("%.3f", Mean)), vjust = -0.8, size = 4) +
  theme_minimal(base_size = 14) +
  scale_y_continuous(
    limits = c(0, 0.6),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = metric_colors) + 
guides(fill = "none") + 
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 15, family = "", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    aspect.ratio = 33/23
  )