library(daprcpp)
library(dplyr)
library(purrr)
library(readr)
library(ldmap)
library(fs)
library(tidyr)
library(stringr)


gf <- snakemake@input[["gwasf"]]
af <- snakemake@input[["annof"]]

torus_cmd <- snakemake@params[["torus_cmd"]]

if (is.null(af)) {
    af <- tempfile()
    write_tsv(tibble::tibble(SNP = character()), af)
}

od <- snakemake@output[["outputf"]]

torus_ret <- daprcpp:::run_torus_single(gf = gf, af = af, torus_path = torus_cmd)
saveRDS(torus_ret, od)
