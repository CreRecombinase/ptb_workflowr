library(EigenH5)
library(tidyverse)


gwas_f <- snakemake@input[["gwasf"]]

snp_f <- snakemake@input[["snpf"]]

gwas_offset <- read_df_h5(gwas_f, "chrom_offset")

snpdb_offset <- read_df_h5(snp_f, "chrom_offset")
