hic_file <- snakemake@input[["hic"]]
    regid <- snakemake@params[["region_id"]]


#    #source("renv/activate.R")
    library(stringr)
    library(readr)
    library(EigenH5)
    library(ldmap)
    library(tidyr)
    library(vroom)
    library(purrr)
    library(dplyr)

    hic_anno_dff <- function(hic_file, ldmrid) {
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
        ldmr <- ldetect_EUR[ldmrid]
        read_tsv(hic_file, col_names = names(cold$cols), col_types = cold, skip = 1L) %>%
            filter(bait_chr != "chrY", otherEnd_chr != "chrY") %>%
            compact_ldmap_region(chrom = "bait_chr", start = "bait_start", end = "bait_end", ldmap_region = "bait") %>%
            compact_ldmap_region(chrom = "otherEnd_chr", start = "otherEnd_start", end = "otherEnd_end", ldmap_region = "target") %>%
            filter(!is.na(region_in_region(bait,ldmr,TRUE)), !is.na(region_in_region(target,ldmr,TRUE)))
    }

    stopifnot(!is.null(regid))

    hic_result_df <- hic_anno_dff(hic_file,as.integer(regid))

    saveRDS(hic_result_df, snakemake@output[["hic"]])
