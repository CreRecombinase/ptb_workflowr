library(EigenH5)
library(readr)
library(ldmap)
library(tidyverse)

mc <- readr::cols(
  "bin" = col_integer(),
  "chrom" = col_character(),
  "pos" = col_double(),
  "chromEnd" = col_integer(),
  "name" = col_character(),
  "score" = col_integer(),
  "strand" = col_factor(levels = c("+", "-")),
  "refNCBI" = col_character(),
  "ref" = col_character(),
  "observed" = col_character(),
  "molType" = col_factor(levels = c(
    "unknown",
    "genomic",
    "cDNA"
  )),
  "class" = col_factor(levels = c(
    "single",
    "in-del",
    "microsatellite",
    "named",
    "mnp",
    "insertion",
    "deletion"
  )),
  "valid" = col_skip(),
  "avHet" = col_skip(),
  "avHetSE" = col_skip(),
  "func" = col_skip(),
  "locType" = col_skip(),
  "weight" = col_skip(),
  "exceptions" = col_skip(),
  "submitterCount" = col_skip(),
  "submitters" = col_skip(),
  "alleleFreqCount" = col_skip(),
  "alleles" = col_skip(),
  "alleleNs" = col_skip(),
  "alleleFreqs" = col_skip(),
  "bitfields" = col_skip()
)

callback_fun <- function(df, filename, datapath, ...) {
  write_df_h5(
    df = dplyr::arrange(
      dplyr::mutate(dplyr::filter(df, class == "single"),
        chrom = fast_str2int(chrom, offset = 3),
        refNCBI = fast_str2ascii(refNCBI),
        ref = fast_str2ascii(ref),
        alt = fast_str2ascii(observed, offset = 2),
        snp_struct = new_ldmap_snp(
          chrom,
          pos,
          ref,
          alt
        ),
        rsid = fast_str2int(name, offset = 2, prefix = "rs")
      ),
      rank.ldmap_snp(snp_struct), ref, alt
    ),
    filename = filename, datapath = datapath, ... = ...
  )
}


input_f <- snakemake@input[["input"]]
output_f <- snakemake@output[["output"]]

file.remove(output_f[file.exists(output_f)])
stopifnot(
  !is.null(input_f),
  !is.null(output_f),
  file.exists(input_f),
  !file.exists(output_f)
)

delim2h5(input_f,
  output_file = output_f,
  h5_args = list(datapath = "snp"),
  delim = "\t",
  col_names = names(mc$cols),
  col_types = mc,
  callback_fun = callback_fun,
  progress = TRUE,
  chunk_size = 300000
)

chrom_vec <- read_vector_h5v(output_f, "snp/chrom", i = integer())
gwas_offset <- rle2offset(chrom_vec) %>%
  rename(chrom = value)
write_df_h5(gwas_offset, output_f, "chrom_offset")
