# install.packages("BiocManager")
BiocManager::install("snpStats")
install.packages("remotes")
remotes::install_github("josefin-werme/LAVA")
library(LAVA)
library(tidyverse)


df <- gwas.ded |>
  select(SNP, A1, A2, NMISS, STAT) |>
  write.table("lava_ded.txt", sep = " ", quote = FALSE, row.names = FALSE)

input = process.input(input.info.file="lava/input_info.txt", # input info file
                      sample.overlap.file=NULL,
                      ref.prefix="lava/g1000_eas/g1000_eas",  # reference genotype data prefix
                      phenos=c("IBS", "DEP", "DED")) # subset of phenotypes listed in the input info file that we want to process

loci = read.loci("lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"); n.loc = nrow(loci)
univ.p.thresh = .05

print(paste("Starting LAVA analysis for",n.loc,"loci"))
progress = ceiling(quantile(1:n.loc, seq(.05,1,.05))) 

# 9, 13, 15, 16, 26, 30-32, 34, 36, 38, 47, 51, 52, 54, 61, 76
# 104, 111, 114, 122, 132-134, 139, 144, 158, 161, 163, 165, 167, 171, 174, 178, 179, 181, 184, 186, 190, 193, 195, 196
# 206, 211, 224, 227, 230, 250, 251, 253, 270, 275, 276, 281, 286, 294, 296, 297
# 313, 316, 328, 333, 334, 338, 354, 359, 360, 362, 374, 385, 298
# 403, 406, 411, 418, 424, 430-432, 435, 436, 444, 449, 450, 455, 461, 471, 478, 480, 484, 499
# 504, 508, 511, 541, 546, 548, 554, 562, 570, 571, 574, 578, 588-591
# 603, 617, 619, 620, 634, 637, 641, 651, 673, 683, 696
# 703, 704, 711, 719, 720, 725, 727, 731-734, 736, 767, 778, 780, 781, 785, 799
# 804, 806, 807, 809, 810, 834, 839, 842, 844, 859, 863, 868, 888
# 907, 912, 913, 917-920, 924, 926, 928, 931, 934, 936, 939, 940, 943, 958, 963, 973, 983, 987, 999
# 1003, 1015, 1033, 1036, 1038, 1043, 1049, 1067, 1074, 1079, 1080, 1083, 1093, 1095, 1099
# 1103, 1115, 1121, 1128, 1133, 1134, 1143, 1146, 1171, 1172, 1182, 1194
# 1214, 1222, 1224, 1225, 1230, 1245, 1251, 1263, 1264, 1266, 1276, 1283
# 1304, 1322, 1325-1327, 1333, 1336, 1351, 1353, 1367, 1376, 1380, 1382, 1395, 1396
# 1401, 1403, 1422, 1424, 1427, 1428, 1438, 1444, 1451, 1453, 1465, 1468, 1484, 1494, 1495
# 1503, 1507, 1511, 1541, 1550, 1552, 1554, 1556, 1565, 1571, 1582, 1588, 1594, 1596-1598
# 1600-1602, 1610, 1614, 1620, 1623, 1630, 1632, 1642, 1646, 1656, 1658, 1659, 1677, 1683, 1692, 1696, 1699
# 1710, 1717, 1718, 1730, 1733, 1735, 1742, 1744, 1745, 1753, 1760, 1765, 1782, 1783, 1790
# 1804, 1816, 1822, 1833, 1835, 1844, 1846, 1848, 1852, 1855, 1856, 1871, 1874, 1875, 1885, 1888, 1894. 1899
# 1918, 1920, 1933, 1937, 1943, 1955, 1956, 1965, 1968, 1986, 1994
# 2003, 2005, 2006, 2009, 2011, 2018-2020, 2025, 2032, 2059, 2062, 2063, 2070, 2078, 2092, 2096
# 2100, 2108, 2109, , 2120, 2124, 2126, 2135, 2136, 2138, 2139, 2143, 2145, 2152, 2155, 2158, 2166, 2169, 2178, 2179, 2181, 2183, 2185-2187, 2196
# 2206, 2207, 2212, 2215, *2224, 2230, 2241, 2244-2248, 2253, 2256, 2257, *2258, 2259, 2265, 2277, 2286, 2287, 2293, 2298
# 2316, 2318, 2327, 2336-2339, 2342, 2347, 2352, 2362, 2371, 2376, 2378, 2390
# 2403-2405, 2408, 2411, 2426, 2439-2441, 2443, 2447, 2449, 2453, 2454, 2463, 2470, 2471, 2473, 2476, 2478, 2482, 2490, 2494

u=b=list()
for (i in 2358:n.loc) {
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
out.fname <- "lava2495"
write.table(do.call(rbind,u), paste0(out.fname,".univ.lava"), row.names=F,quote=F,col.names=T)
write.table(do.call(rbind,b), paste0(out.fname,".bivar.lava"), row.names=F,quote=F,col.names=T)

print(paste0("Done! Analysis output written to ",out.fname,".*.lava"))

