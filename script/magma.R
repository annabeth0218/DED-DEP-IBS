library(tidyverse)
library(data.table)
library(readr)
library(ggvenn)

# for ibs-1g r^2 6e-1
magma_ibs_all <- read.delim("~/R_Projects/ded-dep-ibs/fuma/r2_6e-1/ibs-1g_524698/FUMA_job524698/magma.genes.out") 
magma_ibs <- magma_ibs_all |> filter(P < 0.05)
nrow(magma_ibs) #15763 -> 930 .05 -> 202 .01 -> 110 .005 -> 24 .001
# for ded-1g r^2 6e-1
magma_ded_all <- read.delim("~/R_Projects/ded-dep-ibs/fuma/r2_6e-1/ded-1g_524703/FUMA_job524703/magma.genes.out")
magma_ded <- magma_ded_all |> filter(P < 0.05)
nrow(magma_ded) #15763 -> 1120 .05 -> 273 .01 -> 141 .005 -> 36 .001
# for dep-1g r^2 6e-1
magma_dep_all <- read.delim("~/R_Projects/ded-dep-ibs/fuma/r2_6e-1/dep-1g_524690/FUMA_job524690/magma.genes.out")
magma_dep <- magma_dep_all |> filter(P < 0.05)
nrow(magma_dep) #15759 -> 485 .05 -> 98 .01 -> 50 .005 -> 4 .001

magma.ibsded <- inner_join(magma_ibs, magma_ded, by = "GENE") |>
  select(GENE, SYMBOL = SYMBOL.x, CHR = CHR.x, START = START.x,
         STOP = STOP.x, Nsnps_ibs = NSNPS.x,
         Nsnps_ded = NSNPS.y, P_ibs = P.x, P_ded = P.y)

magma.ibsdep <- inner_join(magma_ibs, magma_dep, by = "GENE") |>
  select(GENE, SYMBOL = SYMBOL.x, CHR = CHR.x, START = START.x,
         STOP = STOP.x, Nsnps_ibs = NSNPS.x,
         Nsnps_dep = NSNPS.y, P_ibs = P.x, P_dep = P.y)

magma.deddep <- inner_join(magma_ded, magma_dep, by = "GENE") |>
  select(GENE, SYMBOL = SYMBOL.x, CHR = CHR.x, START = START.x,
         STOP = STOP.x, Nsnps_ded = NSNPS.x,
         Nsnps_dep = NSNPS.y, P_ded = P.x, P_dep = P.y)

write.table(magma.deddep, "magma_deddep.txt", sep = "\t", row.names = FALSE, quote = FALSE)


magma.genes_filtered <- magma_ibs |>
  filter(P < 0.01)
nrow(magma.genes_filtered)

test <- filter(magma_dep, P < .01)
nrow(test)
test <- filter(test, P < .005)
nrow(test)
test <- filter(test, P < .001)
nrow(test)

ov <- semi_join(magma.genes_filtered, genes.ded.1g, by = c("GENE" = "ensg"))

# MUC16 ENSG00000181143
magma_ded |>
  filter(GENE == "ENSG00000181143") #found
magma_dep |>
  filter(GENE == "ENSG00000181143") #no, P 0.055965
magma_ibs |>
  filter(GENE == "ENSG00000181143") #no, P 0.84454

# TENM2 ENSG00000145934
magma_ded |>
  filter(GENE == "ENSG00000145934") #found
magma_dep |>
  filter(GENE == "ENSG00000145934") #no, 0.83197
magma_ibs |>
  filter(GENE == "ENSG00000145934") #no, 0.82849

# venn
hvh <- list(
  # DED = magma_ded$GENE,
  MAGMA = ov,
  STRING = NodeList$Node
)

# color: "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"
ggvenn(
  hvh, 
  fill_color = c("#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 5
)

same <- intersect(same, magma_ibs$SYMBOL)
write.table(same, 'magma_genes_intersect.txt', quote = FALSE)

g <- 'PAPL'
# NSNPs GPER1 5, PAPL/ACP7 20
filter(magma_ibs, SYMBOL == g) # GPER1 .049 PAPL .045
filter(magma_ded, SYMBOL == g) # GPER1 .048 PAPL .022
filter(magma_dep, SYMBOL == g) # GPER1 .040 PAPL .019

# overlapping genes
magma_ov.dd <- merge(magma_ded, magma_dep[,-c(2,3,4,10)], by = "GENE") #43
magma_ov.id <- merge(magma_ibs, magma_ded[,-c(2,3,4,10)], by = "GENE") #74
magma_ov.ip <- merge(magma_ibs, magma_dep[,-c(2,3,4,10)], by = "GENE") #28
ov <- unique(c(magma_ov.dd$GENE, magma_ov.id$GENE, magma_ov.ip$GENE))
top2 <- unique(magma_ded$GENE[magma_ded$GENE %in% magma_dep$GENE[magma_dep$GENE %in% magma_ibs$GENE]])

ov.string <- unique(ov[ov %in% NodeList$Node]) #18
write.table(ov.string, 'ov_genes.txt', quote = FALSE, row.names = FALSE)


# for STRING
ibs <- data.table(magma_ibs_all)#[order(P)][1:100]
ded <- data.table(magma_ded_all)#[order(P)][1:100]
dep <- data.table(magma_dep_all)#[order(P)][1:100]
ibs[, source := "ibs"]
ded[, source := "ded"]
dep[, source := "dep"]
top.genes.merge <- rbind(ibs, ded, dep)
top.genes.merge <- mutate(top.genes.merge, logP = log(P))
string <- top.genes.merge[order(P), .SD[1], by = GENE]
out <- string298 |> select(GENE)
write.table(out, 'StringList.txt', quote = FALSE, row.names = FALSE)

string298 <- string
string5e2 <- string
stringall <- string

# string298
str.map <- read_tsv('string/top298/string_mapping.tsv')
View(str.map)
df <- merge(string298, str.map, by.x = "GENE", by.y = "queryItem")
string298 <- df[,c(14, 1, 15, 10, 11, 9, 12, 7, 5, 2, 3, 4, 16)]
colnames(string298) <- c("stringId", "ENSG", "stringName", "Symbol", "source", 
                         "P", "logP", "N", "NSNPS", "CHR", "Start", "Stop", "annot")
write_csv(string298, 'full_list.csv')

ref_ibs <- read_excel("ref.xlsx", sheet = "IBS_GWAS", skip = 1)
ref_ibs$our_P <- magma_ibs_all$P[match(ref_ibs$Gene, magma_ibs_all$SYMBOL)]
write.table(ref_ibs, 'ref_ibs.txt', quote = FALSE, row.names = FALSE)