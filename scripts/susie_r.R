library(ldmap)
library(EigenH5)
library(susieR)
library(dplyr)
library(purrr)

ldmrfun <- function(snakemake){
  ld_id <- snakemake@params[["region"]] %||% snakemake@params[["region_id"]]
  if(is.null(ld_id)){
    ichrom <- snakemake@params[["chrom"]]
    istart <- snakemake@params[["start"]]
    iend <- snakemake@params[["end"]]
    stopifnot(!is.null(ichrom),
              !is.null(istart),
              !is.null(iend),
              !is.na(as.integer(ichrom)),
              !is.na(as.integer(istart)),
              !is.na(as.integer(iend)))
    ldmr <- new_ldmap_region(as.integer(ichrom),
                             as.integer(istart),
                             as.integer(iend))
  }else{
    stopifnot(!is.null(ld_id),
              length(ld_id)==1)
    ld_id <- as.integer(ld_id)
    ldmr <- ldetect_EUR[ld_id]
  }
  return(ldmr)
}
ldmr <- ldmrfun(snakemake)



ld_regions <- readRDS(snakemake@input[["ldgf"]]) %>%
  mutate(ldetect_region=ldetect_EUR[region_id]) %>%
  filter(ldetect_region==ldmr)

p <- dim_h5(snakemake@input[["inputf"]], "snp/snp_struct")
gwas_df <- read_df_h5(snakemake@input[["inputf"]],
                      datapath = "snp",
                      subcols = c("snp_struct", "beta", "se", "N"),
                      offset = ld_regions$offset,
                      datasize = ld_regions$datasize
                      ) %>% rename(gwas_snp=snp_struct) %>%
  mutate(snp_pos=clear_alleles(gwas_snp)) %>% filter(!is.na(beta))
ld_df <- read_df_h5(snakemake@input[["ldf"]],"snp") %>%
  rename(ld_snp=snp_struct) %>%
  mutate(ld_id=1:n(),snp_pos=clear_alleles(ld_snp))

run_df <- inner_join(ld_df,gwas_df,by="snp_pos") %>% mutate(am=allele_match(gwas_snp,ld_snp),am=if_else(is.na(am),factor("perfect_match",levels=levels(am)),am),beta=if_else(am=="reverse_match",-beta,beta)) %>% group_by(ld_id)  %>% filter(N==max(N)) %>% slice(1) %>% ungroup() %>% arrange(ld_id)

R <- read_matrix_h5v(snakemake@input[["ldf"]],
                     "R",
                     run_df$ld_id,
                     run_df$ld_id)

run_df <- readRDS(snakemake@input[["priorf"]]) %>%
  select(snp_pos, prior) %>%
  inner_join(run_df)
stopifnot(nrow(run_df)==nrow(R))
h <- 0.1
h_p <- h / p
prior_v <- h_p * nrow(gwas_df)
sres <- susie_suff_stat(
  bhat = run_df$beta,
  shat = run_df$se,
  R = R,
  n = max(run_df$N),
  scaled_prior_variance = prior_v,
  prior_weights = run_df$prior,
  L = 1,
  null_weight = NULL
)
mutate(run_df,pip=susie_get_pip(sres)) %>% 
  saveRDS(snakemake@output[["output_df"]])
saveRDS(sres, snakemake@output[["outputf"]])
