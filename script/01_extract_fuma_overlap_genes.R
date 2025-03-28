library(readxl)
library(dplyr)

# fuma output genes 
file_path <- "fuma/fuma_genes.xlsx"
fuma_genes.ibs <- read_excel(file_path, sheet = "ibs_genes")
fuma_genes.ded <- read_excel(file_path, sheet = "ded_genes")
fuma_genes.dep <- read_excel(file_path, sheet = "dep_genes")

df <- fuma_genes.dep |> 
   filter(posMapSNPs != 0) |> # 5
   filter(eqtlMapSNPs != 0) |> # 6
   filter(ciMap == 'Yes') # 19
  
# deddep : 0
# ibsded : 0
# ibsdep : 10

df <- read_excel('table/Supp_tables.xlsx', sheet = 'eqtl')
write.table(df, 't2eqtl.txt',
            sep = "\t", row.names = FALSE, quote = FALSE)