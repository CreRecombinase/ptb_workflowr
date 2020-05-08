pipf <- snakemake@input[["pipf"]]
    gwasf <- snakemake@input[["gwasf"]]
    ldgf <- snakemake@input[["ldgf"]]
    regid <- snakemake@params[["region_id"]]
    priorf <- snakemake@input[["priorf"]]


#    #source("renv/activate.R")
    library(stringr)
    library(EigenH5)
    library(ldmap)
    library(readr)
    library(tidyr)
    library(vroom)
    library(purrr)
    library(dplyr)



    merge_gwas_susie <- function(gwas_h5, ldetect_rds, susie_rds, ldmrid, prior_rds) {
        stopifnot(
            !is.null(gwas_h5),
            !is.null(ldetect_rds),
            !is.null(susie_rds),
            all(file.exists(c(
                gwas_h5,
                ldetect_rds,
                susie_rds
            )))
        )

        susie_res <- readRDS(susie_rds)
        ldr <- readRDS(ldetect_rds) %>% slice(ldmrid)
        gwas_df <- read_df_h5(gwas_h5,
            "snp",
            subcols = c("snp_struct", "pval"),
            offset = as.integer(ldr$offset),
            datasize = as.integer(ldr$datasize)
        )
        prior_df <- readRDS(prior_rds)
        if (typeof(prior_df$SNP) == "character") {
            prior_df <- prior_df %>%
                separate(SNP, c("chrom", "pos", "ref", "alt"), convert = TRUE) %>%
                dplyr::mutate(chrom = as.integer(stringr::str_replace(chrom, "chr", ""))) %>%
                compact_snp_struct(snp_struct = "SNP")
        }
        ret_match <- dplyr::bind_cols(gwas_df, join_snp(
            gwas_df$snp_struct,
            prior_df$SNP
        )) %>%
            filter(!is.na(index)) %>%
            dplyr::select(snp_struct = match, index, pval)
        ret_match$prior <- prior_df$prior[ret_match$index]

        stopifnot(
            all(snp_in_region(range(ret_match$snp_struct), ldetect_EUR) == ldmrid),
            length(susie_res$alpha) == nrow(ret_match)
        )
        ret_match <- dplyr::select(ret_match, -index) %>%
            dplyr::mutate(pip = susie_res$pip, CS = FALSE)
        if(length(susie_res$sets$cs)!=1)
          return(ret_match)
        stopifnot(length(susie_res$sets$cs) == 1)
        ret_match$CS[susie_res$sets$cs[[1]]] <- TRUE
        return(ret_match)
    }




    stopifnot(!is.null(regid))

    gwas_result_df <- merge_gwas_susie(gwasf, ldgf, pipf, as.integer(regid), priorf)

    saveRDS(gwas_result_df, snakemake@output[["gwas"]])
