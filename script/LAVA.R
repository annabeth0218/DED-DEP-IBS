# install.packages("BiocManager")
BiocManager::install("snpStats")
install.packages("remotes")
remotes::install_github("josefin-werme/LAVA")
library(LAVA)

# input
df <- gwas.dep |> 
  select(SNP, A1, A2, NMISS, STAT, P) |>
  write.table("lava/input/ldsc_ibs.txt", sep = " ", quote = FALSE, row.names = FALSE)

input = process.input(input.info.file="lava/input_info.txt", # input info file
                      sample.overlap.file=NULL,
                      ref.prefix="lava/g1000_eas/g1000_eas",  # reference genotype data prefix
                      phenos=c("IBS", "DEP", "DED")) # subset of phenotypes listed in the input info file that we want to process

# prep loci file
locfile <- read_table('lava/RiskLoci_fuma.txt')
locfile$Nsnps <- 0
for(i in 1:n.loc){
  
  chr <- locfile$CHR[i]  
  start <- locfile$START[i]
  stop <- locfile$STOP[i]
  snps <- c()
  
  df <- gwas.ded |> filter(str_detect(SNP, "^rs")) |>
    filter(CHR_19 == chr & BP_19 >= start & BP_19 <= stop)
  snps <- c(df$SNP, snps)
  df <- gwas.dep |> filter(str_detect(SNP, "^rs")) |>
    filter(CHR_19 == chr & BP_19 >= start & BP_19 <= stop)
  snps <- c(df$SNP, snps)
  df <- gwas.ibs |> filter(str_detect(SNP, "^rs")) |>
    filter(CHR_19 == chr & BP_19 >= start & BP_19 <= stop)
  snps <- c(df$SNP, snps)
  locfile$SNPS[i] <- paste(unique(snps), collapse = ";")
  locfile$Nsnps[i] <- length(unique(snps))
}
write.table(locfile, 'lava/RiskLoci_fuma.txt', sep = '\t', quote = F, row.names = F)

# read loci
loci = read.loci("lava/RiskLoci_fuma.locfile")
n.loc = nrow(loci)
univ.p.thresh = .05

print(paste("Starting LAVA analysis for",n.loc,"loci"))
progress = ceiling(quantile(1:n.loc, seq(.05,1,.05))) 

u=b=list()
for (i in 1:n.loc) {
  if (i %in% progress) print(paste("..",names(progress[which(progress==i)])))     # (printing progress)
  locus = process.locus(loci[i,], input)                                          # process locus
  
  # It is possible that the locus cannot be defined for various reasons (e.g. too few SNPs), so the !is.null(locus) check is necessary before calling the analysis functions.
  if (!is.null(locus)) {
    # extract some general locus info for the output
    loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)
    
    # run the univariate and bivariate tests
    loc.out = run.univ.bivar(locus, univ.thresh = univ.p.thresh)
    u[[i]] = cbind(loc.info, loc.out$univ)
    if(!is.null(loc.out$bivar)) b[[i]] = cbind(loc.info, loc.out$bivar)
  }
}

# save the output
out.fname <- "lava_fumaloci"
write.table(do.call(rbind,u), paste0(out.fname,".univ.lava"), row.names=F,quote=F,col.names=T)
write.table(do.call(rbind,b), paste0(out.fname,".bivar.lava"), row.names=F,quote=F,col.names=T)

print(paste0("Done! Analysis output written to ",out.fname,".*.lava"))

lava2495.b <- read.csv("lava/lava2495.bivar.lava", sep="")
lava2495.u <- read.csv("lava/lava2495.univ.lava", sep="")

lava_fuma.b <- read.csv("lava/lava_fumaloci.bivar.lava", sep="")
lava_fuma.u <- read.csv("lava/lava_fumaloci.univ.lava", sep="")

lava.eda <- filter(lava2495.u, h2.latent >= 1e-2)
# locus 20 - ibs 6 - 1038 p ~ e-6
# locus 21 - ibs 7 - 1263 p ~ e-5



