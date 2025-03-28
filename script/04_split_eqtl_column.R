# Load required package
library(dplyr)
library(readr)
library(stringr)

# Read the file
file_path <- "t2eqtl.txt"
df <- read_tsv(file_path, col_names = c("ENSG", "symbol", "Database_eQTL_Maps"))

# Function to split the third column based on given rules
split_column <- function(x) {
  if (is.na(x) || x == "NA") {
    return(c(x, x))  # Keep NA values unchanged
  } else if (str_count(x, "/") == 1) {
    return(str_split_fixed(x, "/", 2))
  } else if (str_detect(x, "^GTEx/v")) {
    parts <- str_split(x, "/", simplify = TRUE)
    return(c(paste(parts[1], parts[2], sep="/"), parts[3]))
  } else {
    return(c(x, x))
  }
}

# Apply the function to each row
split_results <- t(apply(df["Database_eQTL_Maps"], 1, function(row) split_column(row)))

# Convert the matrix to a data frame
split_df <- as.data.frame(split_results, stringsAsFactors = FALSE)
colnames(split_df) <- c("Database", "eQTL_Maps")

# Combine with original data
df_final <- cbind(df[, 1:2], split_df)

# Save the processed data
output_path <- "t2eqtl_processed.txt"
write.table(df_final, output_path,
            sep = "\t", row.names = FALSE, quote = FALSE)


eqtl_ded.all <- read_table("fuma/fuma_ded/eqtl.txt")
eqtl_dep.all <- read_table("fuma/fuma_dep/eqtl.txt")
eqtl_ibs.all <- read_table("fuma/fuma_ibs/eqtl.txt")

dt <- eqtl_ibs.all |>
  filter(symbol == 'EYS') |>
  filter(tissue != 'EyeGEx' & tissue != 'CRBL')
nrow(dt)