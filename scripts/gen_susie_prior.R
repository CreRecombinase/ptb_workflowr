library(ldmap)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(EigenH5)

gwasf <- snakemake@input[["gwasf"]]
ld_id <- snakemake@params[["region"]]
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
            length(ld_id) == 1)
  ld_id <- as.integer(ld_id)
  ldmr <- ldetect_EUR[ld_id]
}

mvpf <- snakemake@input[["mvp"]]
annof <- snakemake@input[["annof"]]
outf <- snakemake@output[["outf"]]

stopifnot(!is.null(outf),
          !is.null(mvpf),
          !is.null(annof),
          file.exists(annof),
          file.exists(mvpf))

mvd <- readRDS(mvpf) %>%
  unnest(data) %>%
  select(term, estimate)
gwas_df <- ldmap::read_snp_region_h5(gwasf, ldmr, subcols = c("snp_struct")) %>%
  rename(snp_pos = snp_struct) %>%
  mutate(snp_pos=clear_alleles(snp_pos))

annodf <- left_join(gwas_df,
                    read_delim(annof, delim = " ") %>%
                    rename(snp_pos=SNP) %>% 
                    mutate(snp_pos = clear_alleles(ldmap:::parse_ldmap_SNP(snp_pos))))





intercept <- filter(mvd, term == "Intercept") %>%
  pull(estimate)

pivot_longer(annodf,
    cols = c(-snp_pos),
    names_to = c("term"),
    values_to = c("anno")) %>%
    mutate(term = stringr::str_remove(term, "_d$"),
           anno = if_else(is.na(anno), 0, anno)) %>%
    inner_join(mvd) %>%
    group_by(snp_pos) %>%
    summarise(prior = 1 / (1 + exp(- (intercept + sum(anno * estimate))))) %>%
    saveRDS(outf)
