library(tidyverse)

# prep input file
df <- gwas.ibs |> select(CHR, BP) |>
  mutate(V1 = paste0('chr', CHR), V2 = BP) |> select(V1, BP, V2)
write.table(df, 'lava/liftover_ibs.txt', quote = F, row.names = F, sep = ' ', col.names = F)

# read in lifted file
lift <- read_table('lava/hglft_ibs.bed', col_names = c('chr_19_chr', 'bp_19', 'bp')) |>
  mutate(chr_19 = as.numeric(gsub("chr", "", chr_19_chr))) |>
  select(chr_19, bp_19)
lift$bp_19 <- as.numeric(lift$bp_19)

# read in deleted file
df <- read_table('lava/hglft_del_ibs.bed', col_names = F)
lift.fail <- df[seq(2, nrow(df), by = 2), ]
colnames(lift.fail) <- c('chr_38_chr', 'bp_38', 'bp')
lift.fail <- lift.fail |>
  mutate(chr_38 = as.numeric(gsub("chr", "", chr_38_chr))) |>
  select(chr_38, bp_38)
lift.fail$bp_38 <- as.numeric(lift.fail$bp_38)

dt <- gwas.ibs
dt$CHR_19 <- NA
dt$BP_19 <- NA

# Ensure sorting for proper mapping
dt <- dt |> arrange(CHR, BP)
lift <- lift |> arrange(chr_19, bp_19)
lift.fail <- lift.fail |> arrange(chr_38, bp_38)

# Add an index column to ensure correct row-wise mapping
dt$fail <- F
lift.fail$fail <- T

# Mark rows that failed to lift
dt <- dt |>
  left_join(lift.fail, by = c("CHR" = "chr_38", "BP" = "bp_38")) |>
  mutate(fail = coalesce(fail.y, fail.x)) |> select(-fail.x, -fail.y)

lift <- lift |> mutate(index = row_number())
dt <- dt |>
  mutate(index = if_else(fail == FALSE, cumsum(fail == FALSE), NA_integer_))

# Assign lifted coordinates only for rows that are NOT in lift.fail
dt <- dt |>
  left_join(lift, by = "index") |>
  mutate(CHR_19 = ifelse(!is.na(index), chr_19, NA_integer_),
          BP_19 = ifelse(!is.na(index), bp_19, NA_integer_)) |>
   select(-fail, -chr_19, -bp_19, -index)

# Check if NA count matches lift.fail
print(sum(is.na(dt$CHR_19))) # 1528
print(nrow(lift.fail)) # 1518

gwas.ibs <- dt

write.table(dt, 'table/gwas_snp/with_A2/one_sample_gwas/hg19_gwas_ibs.txt',
            sep = '\t', quote = F, row.names = F)
