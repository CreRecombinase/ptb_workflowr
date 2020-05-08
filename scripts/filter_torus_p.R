library(dplyr)
library(purrr)
library(readr)
fdrc <- as.numeric(snakemake@params[["fdrc"]] %||% "0.1")
fdrff <- snakemake@input[["fdrf"]]
readRDS(fdrff) %>% 
  filter(fdr <= fdrc)  %>%
  select(region_id) %>%
  write_tsv(snakemake@output[["off"]], col_names = FALSE)
