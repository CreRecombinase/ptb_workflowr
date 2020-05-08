library(readr)
library(dplyr)
library(purrr)
library(forcats)
library(ldmap)

data(ldetect_EUR)

annof <- snakemake@input[["annot_f"]]
index_f <- snakemake@input[["bimf"]]

anno_n <- snakemake@params[["annot"]]
chromlist <- snakemake@params[["chroms"]]
outputf <- snakemake@output[["outputf"]]

region_l <- purrr::map(annof, ~read_bed(.x)$ldmap_region)
if (length(anno_n) > 0) {
  names(region_l) <- paste0(anno_n,"_d")
  bc <- bim_cols(chrom = col_chromosome(prefix_chr=FALSE))
  index_l <- purrr::walk(index_f, function(x) {
    fe <- file.exists(outputf)
    input_b <- read_plink_bim(x, cols = bc)$snp_struct
    snp_in_regions(input_b, region_l) %>%
      rename(SNP = ldmap_snp) %>%
      write_delim(outputf, delim = " ", append = fe)
  })
}else {
  tibble::tibble(SNP = character()) %>% write_delim(outputf, delim = " ")
}
