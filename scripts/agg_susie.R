cat(paste(unlist(snakemake@input),collapse="\n"))
saveRDS(purrr::map(unlist(snakemake@input),readRDS),snakemake@output[["outputf"]])
