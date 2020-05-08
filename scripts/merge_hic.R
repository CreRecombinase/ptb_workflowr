library(dplyr)
library(forcats)
library(purrr)
library(readr)
library(ldmap)
library(EigenH5)

cold <- cols(
    bait_chr = col_factor(paste0("chr", c(as.character(1:22), c("X", "Y")))),
    bait_start = col_double(),
    bait_end = col_double(),
    bait_name = col_character(),
    otherEnd_chr = col_factor(paste0("chr", c(as.character(1:22), c("X", "Y")))),
    otherEnd_start = col_double(),
    otherEnd_end = col_double(),
    otherEnd_name = col_character(),
    N_reads = col_double(),
    score = col_double()
)
input_hic <- read_tsv(snakemake@input[["inputf"]], col_names = names(cold$cols), col_types = cold, skip = 1L) %>%
    filter(bait_chr != "chrY", otherEnd_chr != "chrY") %>%
    mutate(bait_chr = fct_drop(bait_chr), otherEnd_chr = fct_drop(otherEnd_chr))

baitf <- snakemake@output[["bait"]]
targetf <- snakemake@output[["target"]]
bothf <- snakemake@output[["both"]]

bait_ld <- new_ldmap_region(
    input_hic$bait_chr,
    input_hic$bait_start,
    input_hic$bait_end
)

target_ld <- new_ldmap_region(
    input_hic$otherEnd_chr,
    input_hic$otherEnd_start,
    input_hic$otherEnd_end
)

both_ld <- merge_ldmap_regions(bait_ld, target_ld)

ldmap_region_2_data_frame(bait_ld) %>%
    write_tsv(baitf, col_names = FALSE)
ldmap_region_2_data_frame(target_ld) %>%
    write_tsv(targetf, col_names = FALSE)
ldmap_region_2_data_frame(both_ld) %>%
    write_tsv(bothf, col_names = FALSE)
