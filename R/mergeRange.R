#!/usr/bin/env Rscript
# merge the matched ranges
# Created by Zhili
args = commandArgs(trailingOnly=TRUE)
if(length(args) < 2){
    stop("must have two inputs")
}


output = args[1]
window = as.numeric(args[2]) * 1000

outfile = paste0(output, ".region.txt")

require(data.table)
file1 = paste0(output, "_clump1.clumped")
file2 = paste0(output, "_clump2.clumped")

c1 = data.table()
c2 = data.table()
if(file.exists(file1)){
    c1 = fread(file1)
    if(nrow(c1) != 0){
        c1[, side:=1]
    }
}

if(file.exists(file2)){
    c2 = fread(file2)
    if(nrow(c2) != 0){
        c2[, side:=2]
    }
}

dt.clump = rbind(c1, c2)
dt.clump[, pair:=output]
fwrite(dt.clump, file=paste0(output, ".raw.clump"), sep="\t")

if(length(unique(dt.clump$side)) != 2){
    file.create(outfile, mode="w")
    message("One of the file have no clumped results")
    quit()
}


c1.val = c1[order(CHR, BP), c(1, 3, 4, 5)]
c2.val = c2[order(CHR, BP), c(1, 3, 4, 5)]

setnames(c1.val, "SNP", "SNP1")
setnames(c2.val, "SNP", "SNP2")

c1.val[, low_BP:= BP - window]
c1.val[, up_BP:= BP + window]

dt.t2 = c2.val[c1.val, .(CHR, SNP1, SNP2,BP, low_BP, up_BP), on=.(CHR, BP >= low_BP, BP <= up_BP)]

dt.t2.val = dt.t2[!is.na(SNP2)]

dt.range = dt.t2.val[!duplicated(low_BP)]

chrs = unique(dt.range$CHR)
region_chrs = c()
region_starts = c()
region_ends = c()
top_snps1 = c()
top_snps2 = c()
for(chr.cur in chrs){
    message(chr.cur)
    dt.cur = dt.range[CHR==chr.cur]
    n = nrow(dt.cur)

    sel_idx = 1
    up_BP = dt.cur$up_BP[1]
    while(TRUE){
        bSelect = FALSE
        for(idx in sel_idx:n){
            if(dt.cur[idx]$low_BP > up_BP){
                # prevent from select the present block indicated by idx. 
                idx = idx - 1
                start = dt.cur[sel_idx]$low_BP
                end = dt.cur[idx]$up_BP
                region_chrs = c(region_chrs, chr.cur)
                region_starts = c(region_starts, start)
                region_ends = c(region_ends, end)
                sel_idx = idx + 1
                up_BP = dt.cur$up_BP[sel_idx]
                top_snps1 = c(top_snps1, paste(dt.cur[sel_idx:idx]$SNP1, collapse = ","))
                top_snps2 = c(top_snps2, paste(dt.cur[sel_idx:idx]$SNP2, collapse = ","))
                bSelect = TRUE
                break
            }
        }
        if(idx == n && (!bSelect)){
            start = dt.cur[sel_idx]$low_BP
            end = dt.cur[idx]$up_BP
            region_chrs = c(region_chrs, chr.cur)
            region_starts = c(region_starts, start)
            region_ends = c(region_ends, end)
            top_snps1 = c(top_snps1, paste(dt.cur[sel_idx:n]$SNP1, collapse = ","))
            top_snps2 = c(top_snps2, paste(dt.cur[sel_idx:n]$SNP2, collapse = ","))
            break
        }
    }
}

dt.loc = data.table(trait=output, CHR=region_chrs, start=region_starts, end=region_ends, top_snps1=top_snps1, top_snps2=top_snps2)
dt.loc[, region:=paste0(CHR, "--", start, "-", end)]
setcolorder(dt.loc, c("trait", "region"))
fwrite(dt.loc, file=outfile, sep="\t", col.names=F)

