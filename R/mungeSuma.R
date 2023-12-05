#!/usr/bin/env Rscript
# Formating the summary data
# Created by Zhili
library(optparse)
library(yaml)
library(data.table)
library(stringi)
library(R.utils)


option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, help="input summary", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="output summary", metavar="character"),
    make_option(c("-m", "--map"), type="character", default=NULL, help="mapping summary", metavar="character"),
    make_option(c("-s", "--sample"), type="character", default=NULL, help="mapping summary", metavar="character")
    )

opt_parser = OptionParser(option_list=option_list)

opt = parse_args(opt_parser)

if(is.null(opt$input) | is.null(opt$map) | is.null(opt$output)){
    print_help(opt_parser)
    stop("please fill all the arguments")
}


yaml_head = c("snp", "CHR", "BP", "A1", "A2", "beta", "se","af","p","N","MAFthresh")
info = read_yaml(opt$map)
if(!all(names(info) %in% yaml_head)){
    stop("Map configuration has to contain valid headers")
}

snp = info[["snp"]]
CHR = info[["CHR"]]
BP = info[["BP"]]
A1 = info[["A1"]]
A2 = info[["A2"]]
beta = info[["beta"]]
se = info[["se"]]
af = info[["af"]]
p = info[["p"]]
N = info[["N"]]
maf = as.numeric(info[["MAFthresh"]])
if(is.na(maf)){
    stop("MAFthresh has to be float value range 0 to 1")
}

dt.suma = fread(opt$input)
message("Input summary: ", opt$input)

suma_names = colnames(dt.suma)

if(!CHR %in% suma_names){
    stop("Please specify correct CHR column name")
}

if(!BP %in% suma_names){
    stop("Please specify correct BP column name")
}

if(!A1 %in% suma_names){
    stop("Please specify correct A1 column name")
}

if(!A2 %in% suma_names){
    stop("Please specify correct A2 column name")
}

if(!beta %in% suma_names){
    stop("Please specify correct BETA column name")
}else{
    setnames(dt.suma, beta, "BETA")
}

if(!se %in% suma_names){
    stop("Please specify correct BETA column name")
}else{
    setnames(dt.suma, se, "SE")
}

if(!af %in% suma_names){
    stop("Please specify correct AF column name")
}else{
    setnames(dt.suma, af, "A1_FREQ")
}

if(!p %in% suma_names){
    stop("Please specify correct P column name")
}else{
    setnames(dt.suma, p, "P")
}

n = 0
if(is.null(N)){
    n = countLines(opt$sample)
    dt.suma[, N:=n]
}else{
    n = as.numeric(N)
    if(is.na(n)){
        if(!N %in% suma_names){
            stop("Please specify correct N column name, or a number for sample size")
        }else{
            setnames(dt.suma, N, "N")
        }
    }else{
        dt.suma[, N:=n]
    }
}

snp_sep = stri_match_all(snp, regex="\\{.*?\\}")[[1]][,1]

if(any(is.na(snp_sep))){
    if(length(snp_sep) == 1){
        message("Reconized as SNP column")
        message(snp)
        message(snp_sep)
        message(paste0(suma_names))
        if(!snp %in% suma_names){
            stop("SNP name column is not in the summary data")
        }
        setnames(dt.suma, snp, "SNP")
    }else{
        stop("NA in the SNP names")
    }
}else{
    snp_names = stri_replace_all(snp_sep, regex="[\\{, \\}]", "")
    if(any(!snp_names %in% suma_names)){
        stop("Some columns for SNP name do not exist")
    }
    dt.suma[, SNP:=snp]
    for(idx in 1:length(snp_names)){
        message("Replace SNP name template")
        curCol = snp_names[idx]
        dt.suma$SNP = stri_replace_all(dt.suma$SNP, fixed=snp_sep[idx], dt.suma[[curCol]])
        gc()
    }
    #dt.suma[, SNP:=stri_replace_all(SNP, fixed=snp_sep[1], ..curCol)
}

setnames(dt.suma, c(CHR, BP, A1, A2), c("CHR", "BP", "A1", "A2"))

dt.suma[, MAF:=A1_FREQ]
dt.suma[A1_FREQ > 0.5, MAF := 1 - A1_FREQ]

dt.suma.val = dt.suma[MAF >= maf]
#dt.suma.val[, MAF:=NULL]
dt.suma2 = dt.suma.val[complete.cases(dt.suma.val)]

rm(dt.suma)

saveRDS(dt.suma2, file=paste0(opt$output, ".rds"))
