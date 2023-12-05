#!/usr/bin/env Rscript
# coloc main function
# Created by Zhili
library(coloc)
library(optparse)

option_list = list(
    make_option(c("--finemap1"), type="character", default=NULL, help="finemap1", metavar="character"),
    make_option(c("--finemap2"), type="character", default=NULL, help="finemap2", metavar="character"),
    make_option(c("--output"), type="character", default=NULL, help="output prefix", metavar="character")
    )

opt_parser = OptionParser(option_list=option_list)

opt = parse_args(opt_parser)

if(is.null(opt$finemap1) | is.null(opt$finemap2) | is.null(opt$output)){
    print_help(opt_parser)
    stop("please fill all the arguments")
}


message("Read finemap results")
finemap1 = readRDS(opt$finemap1)
finemap2 = readRDS(opt$finemap2)

message("coloc")
res=coloc.susie(finemap1, finemap2)
saveRDS(res, file=paste0(opt$output, ".rds"))
message("Done!")
