# Install and load necessary packages
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")

library(readxl)
library(dplyr)

# Import 3 sheets from the Excel file
file_path <- "fuma/snp2gene.xlsx"
ibs_data <- read_excel(file_path, sheet = "ibs-1g")
ded_data <- read_excel(file_path, sheet = "ded-1g")
dep_data <- read_excel(file_path, sheet = "dep-1g")

ibs_ded_overlap <- inner_join(ibs_data, ded_data, by = "ensg") %>%
  select(ensg, symbol = symbol.x)

# Overlap between ibs and dep
ibs_dep_overlap <- inner_join(ibs_data, dep_data, by = "ensg") %>%
  select(ensg, symbol = symbol.x)

# Overlap between ded and dep
ded_dep_overlap <- inner_join(ded_data, dep_data, by = "ensg") %>%
  select(ensg, symbol = symbol.x)

# Overlap between all three sheets
all_three_overlap <- reduce(list(ibs_data, ded_data, dep_data), inner_join, by = "ensg") %>%
  select(ensg, symbol = symbol.x)

df <- read_excel(file_path, sheet = "Sheet1")
write.table(df, "ibs_dep_overlap.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# Export the results to text files
write.table(ibs_ded_overlap, "ibs_ded_overlap.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ibs_dep_overlap, "ibs_dep_overlap.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ded_dep_overlap, "ded_dep_overlap.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(all_three_overlap, "all_three_overlap.txt", sep = "\t", row.names = FALSE, quote = FALSE)

dep.fuma <- dep_data |> filter(posMapSNPs != 0) |>
  filter(eqtlMapSNPs != 0)
  # filter(ciMap == 'Yes') # 0

write.table(ibs.fuma, 'ibs_fuma.txt', sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ded.fuma, 'ded_fuma.txt', sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dep.fuma, 'dep_fuma.txt', sep = "\t", row.names = FALSE, quote = FALSE)

