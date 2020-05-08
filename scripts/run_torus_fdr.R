cat("starting!\n")
library(daprcpp)
torus_cmdf <- snakemake@params[["torus_cmd"]]
stopifnot(!is.null(torus_cmdf))
saveRDS(torus_fdr(snakemake@input[["gwasf"]],snakemake@input[["annof"]],torus_path=torus_cmdf),snakemake@output[["outputf"]])
