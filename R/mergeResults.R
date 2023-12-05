#!/usr/bin/env Rscript
# merge all the colocs in a region
# Created by Zhili
args = commandArgs(TRUE)
require(data.table)
require(stringi)

message("Start merging results...")
message(paste0(args, collapse=", "))
traits = args[1]
region = args[2]

out.name = paste0(args[1], "___", args[2], "_combined.tsv")

trait1 = stri_split_fixed(traits, "---", simplify=TRUE)[1]
trait2 = stri_split_fixed(traits, "---", simplify=TRUE)[2]
region = paste0("chr", stri_replace_all(region, fixed="--", ":"))

suma = readRDS(args[3])
coloc = readRDS(args[4])
finemap1 = readRDS(args[5])
finemap2 = readRDS(args[6])

res = coloc$summary

used_snp = unique(c(res$hit1, res$hit2))
idx = match(used_snp, suma$SNP)
suma.val = suma[idx]
rm(suma)

res.out = data.table()

for(curIdx in 1:nrow(res)){
    message("Processing ", curIdx)
    curRes = res[curIdx]
    curIdx1 = curRes$idx1
    curIdx2 = curRes$idx2
    set1 = paste(names(finemap1$sets$cs[[paste0("L", curIdx1)]]), collapse=",")
    set2 = paste(names(finemap2$sets$cs[[paste0("L", curIdx2)]]), collapse=",")
    cover1 = finemap1$sets$coverage[curIdx1]
    cover2 = finemap2$sets$coverage[curIdx2]

    purity1 = finemap1$sets$purity[curIdx1, ]
    purity2 = finemap2$sets$purity[curIdx2, ]

    pip1 = finemap1$pip[curRes$hit1]
    pip2 = finemap2$pip[curRes$hit2]

    suma1 = suma.val[SNP==curRes$hit1]
    chrom1 = suma1$CHR
    BP1 = suma1$BP
    A1_1 = suma1$A1
    A2_1 = suma1$A2
    CHR_BP_A2_A1_1 = paste0("chr", chrom1, "_", BP1, "_", A2_1, "_", A1_1)
    A1_FREQ1 = suma1$A1_FREQ
    BETA1 = suma1$BETA
    SE1 = suma1$SE
    P1 = suma1$P
    N1 = suma1$N
    FREQ_BETA_SE_P_N1 = paste(A1_FREQ1, BETA1, SE1, P1, N1, sep=",")

    suma2 = suma.val[SNP==curRes$hit2]
    chrom2 = suma2$CHR
    BP2 = suma2$BP
    A1_2 = suma2$A1
    A2_2 = suma2$A2
    CHR_BP_A2_A1_2 = paste0("chr", chrom2, "_", BP2, "_", A2_2, "_", A1_2)
    A1_FREQ2 = suma2$A1_FREQ2
    BETA2 = suma2$BETA2
    SE2 = suma2$SE2
    P2 = suma2$P2
    N2 = suma2$N2
    FREQ_BETA_SE_P_N2 = paste(A1_FREQ2, BETA2, SE2, P2, N2, sep=",")

    res.out = rbind(res.out, cbind(data.table(trait1=trait1, trait2=trait2, region=region), curRes, data.table(coverage1=cover1, coverage2=cover2, purity1=paste(purity1, collapse=","), purity2=paste(purity2, collapse=","), pip1=pip1, pip2=pip2, CHR_BP_A2_A1_1=CHR_BP_A2_A1_1, CHR_BP_A2_A1_2=CHR_BP_A2_A1_2, FREQ_BETA_SE_P_N1=FREQ_BETA_SE_P_N1, FREQ_BETA_SE_P_N2=FREQ_BETA_SE_P_N2,set1=set1, set2=set2)))
}


fwrite(res.out, file=out.name, sep="\t", na="NA", quote=F)
message("Done")
