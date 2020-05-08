library(EigenH5)
library(ldmap)
library(dplyr)
library(purrr)

data(ldetect_EUR)
rdsf <- unlist(snakemake@input)
stopifnot(!is.null(rdsf), file.exists(rdsf))

ldmrfun <- function(snakemake) {
    ld_id <- snakemake@params[["region"]]
    if (is.null(ld_id)) {
        ichrom <- snakemake@params[["chrom"]]
        istart <- snakemake@params[["start"]]
        iend <- snakemake@params[["end"]]
        stopifnot(
            !is.null(ichrom),
            !is.null(istart),
            !is.null(iend),
            !is.na(as.integer(ichrom)),
            !is.na(as.integer(istart)),
            !is.na(as.integer(iend))
        )
        ldmr <- new_ldmap_region(
            as.integer(ichrom),
            as.integer(istart),
            as.integer(iend)
        )
    } else {
        stopifnot(!is.null(ld_id), length(ld_id) == 1)
        ld_id <- as.integer(ld_id)
        ldmr <- ldetect_EUR[ld_id]
    }
    return(ldmr)
}
ldmr <- ldmrfun(snakemake)


init_fn <- function(reference_files, ldmr) {
    return(reference_files)
}
filter_map_fn <- function(reference_file, ldmr) {
    dplyr::filter(read_df_h5(reference_file, "snp"), snp_struct %overlaps% ldmr)
}
filter_geno_fn <- function(reference_file, ldmr, map) {
    read_matrix_h5v(reference_file, "dosage", j = map$index)
}
write_fn <- function(map, bsx, output_file) {
    write_df_h5(map, output_file, "snp")
    write_matrix_h5(bsx, output_file, "dosage")
    return(output_file)
}
read_map_fn <- function(x) {
    read_df_h5(x, "snp")
}
read_dosage_fn <- function(x, ...) {
    read_matrix_h5v(x, "dosage")
}
srds <- ldmap::subset_rds(
    ldmr = ldmr,
    reference_files = rdsf,
    output_file = snakemake@output[["ldf"]],
    init_fn = init_fn,
    filter_map_fn = filter_map_fn,
    filter_geno_fn = filter_geno_fn,
    write_fn = write_fn
)
print(srds)
retR <- panel_ld(srds,
    LDshrink = FALSE,
    read_map_fn = read_map_fn,
    read_dosage_fn = read_dosage_fn
)
attr(retR, "dimnames") <- NULL
EigenH5::write_matrix_h5(retR, snakemake@output[["ldf"]], datapath = "R")
