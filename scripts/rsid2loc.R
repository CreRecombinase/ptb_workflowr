library(dplyr)
library(purrr)
library(readr)


library(ldmap)


  input_f <- snakemake@input[["input"]]
  output_f <- snakemake@output[["output"]]
  input_ids <- EigenH5::fast_str2int(scan(input_f, what = character()), prefix = "rs")
  input_ids <- input_ids[!is.na(input_ids)]
  BSgenome::snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37::SNPlocs.Hsapiens.dbSNP144.GRCh37,
                     ids = input_ids,
                     ifnotfound = "warn") %>% as_tibble() %>% 
      dplyr::rename(chrom = seqnames, rsid = RefSNP_id) %>%
      dplyr::mutate(chrom = as.integer(chrom),
                    rsid = rsid) %>%
      select(-strand) %>%
      readr::write_tsv(output_f)
