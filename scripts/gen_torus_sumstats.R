library(readr)
library(dplyr)
library(purrr)
library(forcats)
library(ldmap);
library(EigenH5)


data(ldetect_EUR)


sumstat_h5f <- snakemake@input[["inputf"]]
snplist <- snakemake@input[["snplist"]]
chromlist <- snakemake@params[["chroms"]]
outputf <- snakemake@output[["outputf"]]

chrom_df <- read_df_h5(sumstat_h5f, "chrom_offset") %>% 
  dplyr::slice(1:22) %>% 
  dplyr::mutate(offset = as.integer(offset),
                datasize = as.integer(datasize)) %>%
  dplyr::arrange(offset)

bc <- bim_cols(chrom=col_chromosome(prefix_chr=FALSE))
index_l <- purrr::map(snplist, ~read_plink_bim(.x,cols = bc)$snp_struct)

write_snplist_txt <- function(chrom, offset, datasize, snplist_l, ...) {
      stopifnot(inherits(snplist_l, "ldmap_snp"))
      fe <- file.exists(outputf)
      idf <- EigenH5::read_df_h5(
                        filename = sumstat_h5f,
                        datapath = "snp",
                        subcols = c("snp_struct", "beta", "se"),
                        offset = offset,
                        datasize = datasize) %>%
        mutate(ref_snp = snplist_l[snp_in_snp(snp_struct, snplist_l)],
               match_type = allele_match(snp_struct, ref_snp),
               beta = dplyr::if_else(is_reversed(match_type),
                                     -beta,
                                     beta)) %>% 
        filter(!is.na(ref_snp)) %>% 
        dplyr::transmute(SNP = ref_snp,
                         locus = ldmap:::snp_in_region(ref_snp, ldetect_EUR),
                         `z-vals` =  beta / se )
      stopifnot(all(!is.na(idf$locus)))
      write_delim(idf, outputf, delim = " ", append = fe)
    }


mutate(chrom_df, snplist_l = index_l) %>%
  pwalk(write_snplist_txt)
