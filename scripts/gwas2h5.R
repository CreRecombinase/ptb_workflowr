library(dplyr)
library(purrr)
library(readr)
library(EigenH5)
library(readr)
library(ldmap)

input_f <- snakemake@input[["inputf"]]
output_f <- snakemake@output[["outputf"]]
paramf <- snakemake@input[["paramf"]]
stopifnot(!is.null(paramf))
source(paramf)

callback_fun <- function(df, filename, datapath, ...){
  write_df_h5(
    df = dplyr::slice(
                  dplyr::mutate(df,
                                ref = fast_str2ascii(A2),
                                alt = fast_str2ascii(A1),
                                snp_struct =
                                  new_ldmap_snp(chrom, pos, ref, alt),
                                rsid = fast_str2int(rsid, prefix = "rs"),
                                ),
                  rank.ldmap_snp(snp_struct)),
    filename = filename, datapath = datapath, ... = ...)
}

stopifnot(!is.null(input_f),
          !is.null(output_f),
          file.exists(input_f),
          !file.exists(output_f))

delim2h5(input_f,
         output_file = output_f,
         h5_args = list(datapath = "snp"),
         delim = data_delim,
         col_names = names(mc$cols),
         skip = 1L,
         callback_fun = callback_fun,
         col_types = mc,
         progress = TRUE,
         chunk_size = 150000)

chrom_vec <- read_vector_h5v(output_f, "snp/chrom", i = integer())
chrom_df <- rle2offset(as.integer(chrom_vec)) %>%
  dplyr::rename(chrom = value) %>% 
  mutate(offset=as.integer(offset),datasize=as.integer(datasize))
write_df_h5(chrom_df,output_f,"chrom_offset")
