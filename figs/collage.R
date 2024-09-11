library(pdftools)
library(gridExtra)
library(ggplot2)
library(magick) # more versatile with pdf
library(TwoSampleMR)
library(png)

# pdf to grobs
pdf_to_grob <- function(file) {
  img <- image_read_pdf(file, density = 300)
  img <- image_trim(img)  # Trim any excess white space
  g <- rasterGrob(img, interpolate = TRUE)
  return(g)
}

pdf_dir <- 'loo_ibs_dry_dep/LOO_pdf/two'

plot_order <- c(
  'leaf_one_out_expibs_outdep_twoadj.pdf' = 1,
  'leaf_one_out_expdep_outibs_twoadj.pdf' = 2,
  'leaf_one_out_expibs_outdry_twoadj.pdf' = 3,
  'leaf_one_out_expdry_outibs_twoadj.pdf' = 4,
  'leaf_one_out_expdry_outdep_twoadj.pdf' = 5,
  'leaf_one_out_expdep_outdry_twoadj.pdf' = 6,
  'leaf_one_out_expibs_outdep_twocojo.pdf' = 7,
  'leaf_one_out_expdep_outibs_twocojo.pdf' = 8,
  'leaf_one_out_expibs_outdry_twocojo.pdf' = 9,
  'leaf_one_out_expdry_outibs_twocojo.pdf' = 10,
  'leaf_one_out_expdry_outdep_twocojo.pdf' = 11,
  'leaf_one_out_expdep_outdry_twocojo.pdf' = 12,
  'leaf_one_out_expibs_outdep_twosaige.pdf' = 13,
  'leaf_one_out_expdep_outibs_twosaige.pdf' = 14,
  'leaf_one_out_expibs_outdry_twosaige.pdf' = 15,
  'leaf_one_out_expdry_outibs_twosaige.pdf' = 16,
  'leaf_one_out_expibs_outdep_twosaige.pdf' = 17, # false
  'leaf_one_out_expibs_outdry_twosaige.pdf' = 18, # false
  'leaf_one_out_expibs_outdep_twomtag.pdf' = 19,
  'leaf_one_out_expdep_outibs_twomtag.pdf' = 20,
  'leaf_one_out_expibs_outdry_twomtag.pdf' = 21,
  'leaf_one_out_expdry_outibs_twomtag.pdf' = 22,
  'leaf_one_out_expdry_outdep_twomtag.pdf' = 23,
  'leaf_one_out_expdep_outdry_twomtag.pdf' = 24
)

pdf_files <- file.path(pdf_dir, names(sort(plot_order)))
plot_list <- lapply(pdf_files, pdf_to_grob)
plot_list[[17]] <- nullGrob()
plot_list[[18]] <- nullGrob()

n_plots <- length(plot_list)
n_cols <- 6
n_rows <- ceiling(n_plots / n_cols)
combined_plot <- do.call(arrangeGrob, c(plot_list, list(nrow = n_rows, 
                                                        ncol = n_cols,
                                                        widths = unit(rep(1, n_cols), "null"),
                                                        heights = unit(rep(1, n_rows), "null"),
                                                        padding = unit(0.5, "line"))))
ggsave("combined_plots.pdf", combined_plot, width = 30, height = 20, 
       units = "in", device = cairo_pdf, dpi = 300)

# pdf to png
pdf_dir <- 'MR_scatter_plots/exp_DEP_out_IBS_mr.pdf'
pdf_convert(pdf_dir, format = "png", dpi = 300)