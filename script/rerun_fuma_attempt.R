library(tidyverse)
library(readxl)

# prepare magma annot

# --1: prep snploc--
snps <- read.delim("table/gwas_snp/with_A2/liftover/hg19_gwas_ded.txt")
s.out <- snps |> select(SNP, CHR_19, BP_19)
write.table(s.out, 'fuma/rerun/snploc.txt', sep = "\t", row.names = FALSE, quote = FALSE)

# --2: prep geneloc--
genes <- read.delim("fuma/fuma_ded/magma.genes.out")
g.out <- genes |> filter(P < 0.05) |> select(GENE, CHR, START, STOP)
write.table(g.out, 'fuma/rerun/geneloc.txt', sep = "\t", row.names = FALSE, quote = FALSE)
            
# prepare fuma
fuma.out <- snps |> select(SNP, CHR_19, BP_19, A1, A2, P, OR, SE, NMISS)
write.table(fuma.out, 'fuma/rerun/fuma_input_ded.txt', sep = "\t", row.names = FALSE, quote = FALSE)

# prepare HiC
GSE192625_Hi_C_interaction <- read_excel("/Volumes/Anna_256/HiC/GSE192625_Hi-C_interaction.xlsx", 
                                         sheet = "cis_interaction.annotation")

hic.out <- GSE192625_Hi_C_interaction |>
  separate(`chr1:start1:end1`, into = c("chr1", "start1", "end1"), sep = "[:-]") |>
  separate(`chr2:start2:end2`, into = c("chr2", "start2", "end2"), sep = "[:-]") |>
  select(chr1, start1, end1, chr2, start2, end2) |>
  mutate(fdr = 0)

write.table(fuma.out, '/Volumes/Anna_256/HiC/GSE192625.txt', sep = "\t", row.names = FALSE, quote = FALSE)

# plot manhattan
library(qqman)
library("CMplot")

input <- gwas.ibs |> select(SNP, CHR, BP, P)
colnames(input) <- c("SNP", "chr", "pos", "IBS")

CMplot(input, plot.type="m", col=c("grey30","grey60"), 
       LOG10=TRUE, ylim=c(0,7), cex=0.6, chr.labels = 1:22,
       threshold=1e-5,threshold.lty=2, threshold.lwd=1, threshold.col="grey", 
       amplify=TRUE,chr.den.col=NULL, signal.col="#BA8E23", signal.cex=1.5,signal.pch=19,
       file="jpg",file.name="ibs",dpi=300,
       file.output=TRUE,verbose=TRUE,width=14,height=6)

input <- gwas.ded |> select(SNP, CHR, BP, P)
colnames(input) <- c("SNP", "chr", "pos", "IBS")

CMplot(input, plot.type="m", col=c("grey30","grey60"), 
       LOG10=TRUE, ylim=c(0,7), cex=0.6, chr.labels = 1:22,
       threshold=1e-5,threshold.lty=2, threshold.lwd=1, threshold.col="grey", 
       amplify=TRUE,chr.den.col=NULL, signal.col="#950606", signal.cex=1.5,signal.pch=19,
       file="jpg",file.name="ded",dpi=300,
       file.output=TRUE,verbose=TRUE,width=14,height=6)

input <- gwas.dep |> select(SNP, CHR, BP, P)
colnames(input) <- c("SNP", "chr", "pos", "IBS")

CMplot(input, plot.type="m", col=c("grey30","grey60"), 
       LOG10=TRUE, ylim=c(0,7), cex=0.6, chr.labels = 1:22,
       threshold=1e-5,threshold.lty=2, threshold.lwd=1, threshold.col="grey", 
       amplify=TRUE,chr.den.col=NULL, signal.col="#000080", signal.cex=1.5,signal.pch=19,
       file="jpg",file.name="dep",dpi=300,
       file.output=TRUE,verbose=TRUE,width=14,height=6)

input <- gwas.dep |>
  select(SNP, CHR, BP, DEP = P) |>
  inner_join(gwas.ded |> select(SNP, DED = P), by = "SNP") |>
  inner_join(gwas.ibs |> select(SNP, IBS = P), by = "SNP")

CMplot(input,type="p",plot.type="c",r=0.4,col=c("grey30","grey60"), cex=0.4,
         chr.labels=paste("Chr",c(1:22),sep=""), ylim=c(0,7),
         threshold=c(1e-5, 1e-4),threshold.lty=c(1, 2), threshold.lwd=c(1, 1), threshold.col=c("red","blue"),
         amplify=TRUE,chr.den.col="black", signal.col=c("red","blue"), signal.cex=0.8,signal.pch=19,
         bin.size=1e6,outward=FALSE,
         file="jpg",file.name=NULL,dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)

