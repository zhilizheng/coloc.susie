#!/usr/bin/env Rscript
# merge the summary data, keep the common variants only
# Created by Zhili
library(data.table)
library(stringi)

args = commandArgs(trailingOnly=TRUE)

file1 = paste0(args[1], ".rds")
file2 = paste0(args[2], ".rds")

p1 = as.numeric(args[3])
p2 = as.numeric(args[4])

output = args[5]

dt1 = readRDS(file1)
dt2 = readRDS(file2)

idx = match(dt2$SNP, dt1$SNP)
idx2 = which(is.finite(idx))

dt1.com = dt1[idx[idx2]]
dt2.com = dt2[idx2]

message(nrow(dt1.com), " SNPs in common")

if(nrow(dt1.com) != sum(dt1.com$SNP == dt2.com$SNP)){
    stop("Something unexpeced happened after matching the SNPs")
}

bA1A1 = (dt1.com$A1 == dt2.com$A1) & (dt1.com$A2 == dt2.com$A2)
bA1A2 = (dt1.com$A1 == dt2.com$A2) & (dt1.com$A2 == dt2.com$A1)

dt2.com[bA1A2, A1_FREQ := 1 - A1_FREQ]
dt2.com[bA1A2, BETA:= (-BETA)]

dt1.com[, A1_FREQ2:=dt2.com$A1_FREQ]
dt1.com[, MAF2:=dt2.com$MAF]
dt1.com[, BETA2:=dt2.com$BETA]
dt1.com[, SE2:=dt2.com$SE]
dt1.com[, P2:=dt2.com$P]
dt1.com[, N2:=dt2.com$N]

dt1.com.val = dt1.com[bA1A1 | bA1A2]
rm(dt1.com, dt2.com, dt1, dt2)
gc()

# for sex name to X (consistent with genotype)
dt1.com.val[, SNP:=stri_replace_all(SNP, fixed="chr23", "chrX")]
saveRDS(dt1.com.val, file=paste0(output, ".rds"))

# exclude MHC
idx = dt1.com.val[!(CHR=="6" & BP >= 28510120 & BP <= 33480577), which=TRUE]
## suma 1
dt.out = dt1.com.val[idx, .(CHR, BP, SNP, A1, A2, P)][P <= p1]
fwrite(dt.out, file=paste0(output, "_clump1.txt"), sep="\t")
cat(dt.out$SNP, file=paste0(output, "_clump1.snplist"), sep="\n")

## suma2
dt.out = dt1.com.val[idx, .(CHR, BP, SNP, A1, A2, P2)][P2 <= p2]
setnames(dt.out, "P2", "P")
fwrite(dt.out, file=paste0(output, "_clump2.txt"), sep="\t")
cat(dt.out$SNP, file=paste0(output, "_clump2.snplist"), sep="\n")
