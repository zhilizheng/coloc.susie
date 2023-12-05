#!/usr/bin/env Rscript
# prepare for the coloc
# Created by Zhili
library(optparse)
library(yaml)
library(data.table)
library(stringi)

option_list = list(
    make_option(c("--suma"), type="character", default=NULL, help="input summary", metavar="character"),
    make_option(c("--prefix"), type="character", default=NULL, help="output prefix", metavar="character"),
    make_option(c("--region"), type="character", default=NULL, help="string indicate region, CHR:start-end", metavar="character"),
    make_option(c("--bgenInf1"), type="character", default=NULL, help="bgen information for 1, yaml", metavar="character"),
    make_option(c("--bgenInf2"), type="character", default=NULL, help="bgen information for 2, yaml", metavar="character")
    )

opt_parser = OptionParser(option_list=option_list)

opt = parse_args(opt_parser)

if(is.null(opt$prefix) | is.null(opt$region) | is.null(opt$suma) | is.null(opt$bgenInf1) | is.null(opt$bgenInf2)){
    print_help(opt_parser)
    stop("please fill all the arguments")
}


yaml_head = c("bgen", "bgi", "sample")
info1 = read_yaml(opt$bgenInf1)
info2 = read_yaml(opt$bgenInf2)
if(!all(names(info1) %in% yaml_head)){
    stop("bgen information 1 has to contain valid keys")
}

if(!all(names(info2) %in% yaml_head)){
    stop("bgen information 2 has to contain valid keys")
}

prefix = opt$prefix
region = opt$region

message("Processing region text: ", region)
regions = stri_split_fixed(region, "--", simplify=TRUE)
CHR = regions[1]
poses = stri_split_fixed(regions[2], "-", simplify=TRUE)
start = as.numeric(poses[1])
end = as.numeric(poses[2])
message(" CHR: ", CHR, ", start: ", start, ", end: ", end)

#### output the bgen information
bgen1 = stri_replace_all(info1$bgen, fixed="{CHR}", CHR)
bgi1 = stri_replace_all(info1$bgi, fixed="{CHR}", CHR)
sample1 = stri_replace_all(info1$sample, fixed="{CHR}", CHR)
cat(bgen1, file=paste0(prefix, ".bgenGeno1.txt"))
cat(bgi1, file=paste0(prefix, ".bgenIndex1.txt"))
cat(sample1, file=paste0(prefix, ".bgenSample1.txt"))

bgen2 = stri_replace_all(info2$bgen, fixed="{CHR}", CHR)
bgi2 = stri_replace_all(info2$bgi, fixed="{CHR}", CHR)
sample2 = stri_replace_all(info2$sample, fixed="{CHR}", CHR)
cat(bgen2, file=paste0(prefix, ".bgenGeno2.txt"))
cat(bgi2, file=paste0(prefix, ".bgenIndex2.txt"))
cat(sample2, file=paste0(prefix, ".bgenSample2.txt"))

dt.suma = readRDS(opt$suma)

curCHR = CHR
dt.cur = dt.suma[CHR==curCHR & BP >= start & BP <= end]

message("Processing summary data")
dt.pos = dt.cur[, .(SNP, CHR, BP, A2, A1)]
setnames(dt.pos, c("rsid", "chromosome", "position", "allele1", "allele2"))
dt.pos[, chromosome:=paste0("chr", chromosome)]
zfile = paste0(prefix, ".z")
fwrite(dt.pos, file=zfile, sep=" ")

message("Output the sample size")
N1 = median(dt.cur$N)
N2 = median(dt.cur$N2)
cat(N1, file=paste0(prefix, ".N1.txt"))
cat(N2, file=paste0(prefix, ".N2.txt"))

nsnps = nrow(dt.cur)
cat(nsnps, file=paste0(prefix, ".n.txt"))

D1 = list()
D1[["beta"]] = dt.cur$BETA
D1[["varbeta"]] = dt.cur$SE^2
D1[["N"]] = median(dt.cur$N)
D1[["type"]] = "quant"
D1[["MAF"]] = dt.cur$MAF
D1[["snp"]] = dt.cur$SNP
D1[["position"]] = dt.cur$BP
saveRDS(D1, file=paste0(prefix, ".suma1.rds"))

D2 = list()
D2[["beta"]] = dt.cur$BETA2
D2[["varbeta"]] = dt.cur$SE2^2
D2[["N"]] = median(dt.cur$N2)
D2[["type"]] = "quant"
D2[["MAF"]] = dt.cur$MAF2
D2[["snp"]] = dt.cur$SNP
D2[["position"]] = dt.cur$BP
saveRDS(D2, file=paste0(prefix, ".suma2.rds"))

message("colocPrep done")

