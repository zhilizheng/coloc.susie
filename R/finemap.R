#!/usr/bin/env Rscript
# finemapping 
# Created by Zhili
library(coloc)
library(data.table)
library(optparse)


option_list = list(
    make_option(c("--suma"), type="character", default=NULL, help="input summary", metavar="character"),
    make_option(c("--ldz"), type="character", default=NULL, help="input LD", metavar="character"),
    make_option(c("--output"), type="character", default=NULL, help="output prefix", metavar="character")
    )

opt_parser = OptionParser(option_list=option_list)

opt = parse_args(opt_parser)

if(is.null(opt$suma) | is.null(opt$ldz) | is.null(opt$output)){
    print_help(opt_parser)
    stop("please fill all the arguments")
}


D = readRDS(opt$suma)
dt.ld = fread(opt$ldz)
D[["LD"]] = as.matrix(dt.ld)

colnames(D[["LD"]]) = D[["snp"]]
rownames(D[["LD"]]) = D[["snp"]]

rm(dt.ld)
gc()

message("Run susie")
finemap = runsusie(D)
message("Save results")
saveRDS(finemap, file=paste0(opt$output, ".rds"))
message("Done!")
