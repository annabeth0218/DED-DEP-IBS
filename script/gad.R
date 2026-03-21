library(readr)
alzheimer <- read_csv("~/Library/Mobile Documents/com~apple~CloudDocs/lab/25 HK CityU/cleaning/alzheimer/alzheimer.csv")

library(ggplot2)

# 設定存圖資料夾
output_dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/lab/25 HK CityU/cleaning/alzheimer"  # 自訂你的路徑
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 計算 k
k <- ((ncol(alzheimer) - 1) %/% 2)

for (i in seq(3, 2 * k + 1, by = 2)) {
  # 欄位名稱
  name_col <- colnames(alzheimer)[i]
  count_col <- colnames(alzheimer)[i + 1]
  
  # 擷取 study count 欄，計算 frequency
  counts <- alzheimer[[count_col]]
  counts_table <- as.data.frame(table(counts))
  colnames(counts_table) <- c("study_count", "num_patients")
  counts_table$study_count <- as.numeric(as.character(counts_table$study_count))
  
  # 畫 bar plot + 標註高度
  p <- ggplot(counts_table, aes(x = study_count, y = num_patients)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = num_patients), vjust = -0.5, size = 3) +
    labs(title = paste("Study count for", name_col),
         x = "Study count",
         y = "Number of patients") +
    theme_minimal()
  
  # 存檔名：取 '/' 後的部分
  save_name <- sub(".*/", "", name_col)
  save_path <- file.path(output_dir, paste0(save_name, ".png"))
  
  # 儲存圖片
  ggsave(filename = save_path, plot = p, width = 6, height = 4, dpi = 300)
}

library(ggvenn)

G_bypatient <- read_csv("~/Library/Mobile Documents/com~apple~CloudDocs/lab/25 HK CityU/cleaning/G_bypatient.csv")
LTG_bypatient <- read_csv("~/Library/Mobile Documents/com~apple~CloudDocs/lab/25 HK CityU/cleaning/LTG_bypatient.csv")

venn3 <- list(
  AD_txt = alzheimer$alzheimer_all_txt,
  LTG = gsub("-", "", LTG_bypatient$patient_id[LTG_bypatient$ct != 0]),
  G = gsub("-", "", G_bypatient$patient_id[G_bypatient$ct != 0])
  # AD_2011 = gsub("-", "", alzheimer$alzheimer_2011_bmp),
  # AD_2011_G  = gsub("-", "", alzheimer$`alzheimer_2011-Glaucoma_bmp`)
)

# color: "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"
ggvenn(
  venn3, 
  fill_color = c("#0073C2FF", "#CD534CFF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 5, text_size = 3.6
)

venn2 <- list(
  AD_2011 = alzheimer$alzheimer_2011_bmp,
  AD_2011_G = alzheimer$`alzheimer_2011-Glaucoma_bmp`
)

ggvenn(
  venn2, 
  fill_color = c("#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 5
)
