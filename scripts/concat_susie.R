library(purrr)
in_f <- snakemake@input[["input_rds"]]
saveRDS(map_df(in_f,readRDS),snakemake@output[["output_rds"]])
