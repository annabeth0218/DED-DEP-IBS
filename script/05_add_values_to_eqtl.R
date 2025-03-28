# Load required library
library(dplyr)
library(readr)

# Read input data
t2eqtl <- read_excel('table/Supp_tables.xlsx', sheet = '4') |>
  select(ENSG, symbol, Database, eQTL_Maps)
eqtl <- read_tsv("fuma/fuma_ded/eqtl.txt", col_names = TRUE, skip = 0) |> 
  filter(!is.na(RiskIncAllele))

# Rename columns to avoid issues with special characters
colnames(t2eqtl) <- c("ENSG", "symbol", "Database", "eQTL_Maps")
colnames(eqtl) <- c("uniqID", "db", "tissue", "gene", "testedAllele", "p", "signed_stats", "FDR", "RiskIncAllele", 
                    "alignedDirection", "chr", "pos", "symbol_eqtl", "eqtlMapFilt")

# Convert p-values to numeric (if they are not already)
eqtl$p <- as.numeric(eqtl$p)

# Ensure t2eqtl has enough rows
if (nrow(t2eqtl) >= 288) {
  for (i in 2:281) {
    row_x <- t2eqtl[i, ] # Extract current row
    
    # Filter eqtl rows that match criteria
    matched_rows <- eqtl %>%
      filter(gene == row_x$ENSG, 
             db == row_x$Database, 
             tissue == row_x$eQTL_Maps)
    
    # Select row with smallest p-value
    if (nrow(matched_rows) > 0) {
      best_match <- matched_rows %>%
        arrange(p) %>%
        slice(1) # Get the row with smallest p-value
      
      # Append selected row data as new columns
      t2eqtl[i, paste0("matched_", names(best_match))] <- best_match
    }
  }
}

# Save the modified dataset
output_file <- "t2eqtl_02132354.txt"
t2eqtl <- t2eqtl |> filter(!is.na(matched_uniqID))
write_tsv(t2eqtl, output_file)

# Display message with output file path
cat("Updated file saved at:", output_file)



g <- 'C10orf115'

df <- eqtl |>
  filter(symbol == g) |>
  filter(db == 'PsychENCODE') |>
  filter(tissue == 'PsychENCODE_eQTLs')
View(df)

df <- eqtl |>
  filter(symbol == g) |>
  filter(db == 'eQTLGen') |>
  filter(tissue == 'eQTLGen_cis_eQTLs')
View(df)

df <- eqtl |>
  filter(db == 'CMC') |>
  filter(!is.na(p))
View(df)

