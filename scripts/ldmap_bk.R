source("renv/activate.R")
library(bigsnpr)
rds <- snp_readBed(snakemake@input[["bedlist"]], backingfile = snakemake@params[["rdsp"]])
