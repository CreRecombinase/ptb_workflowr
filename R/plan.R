ld_df <-read_tsv(ldetect_file)
t_ld_df <- mutate(ld_df) %>% select(-region_id)
plan <- drake_plan(
  snp_df = target(
    read_df_reg("EUR",chrom,start,stop),
                  transform = map(.data=!!t_ld_df)),
  mapd_df = target(read_map(chrom,start,stop),
                   transform = map(.data = !!t_ld_df)),
  gwas_df = target(
    snp_reg(t_chrom = chrom,
            t_start = start,t_stop = stop,beta_v = "beta",se_v = "se",N_v="N"),
    transform = map(.data=!!ld_df)
  ),
  input_d =target(merge_df(snp_df = snp_df,gwas_df = gwas_df) ,
                  transform = map(snp_df,gwas_df)),
    mmap_df = target(merge_map(input_d,mapd_df),
                     transform = map(input_d,mapd_df)),
  est_df =target(local_rssp_est(mmap_df,"EUR"),transform = map(mmap_df)),
  full_est_df = target(bind_rows(est_df),transform = combine(est_df)),
  trace=T,
)
