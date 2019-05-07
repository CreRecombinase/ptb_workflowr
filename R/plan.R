ld_df <-read_tsv(ldetect_file)
t_ld_df <- mutate(ld_df) %>% select(-region_id)
ld_df_l <- map(parallel::splitIndices(nrow(ld_df),50),~slice(ld_df,.x))
ind_i <-seq_along(ld_df_l)
plan <- drake_plan(
  snp_df = target(
    map_df_reg("EUR",ld_df_l[[i]]),
    transform = map(i = !!ind_i)),
  mapd_df = target(map_read_map(ld_df_l[[i]]),
                   transform = map(i = !!ind_i)),
  gwas_df = target(
    map_snp_reg(ld_df_l[[i]], beta_v = "beta", se_v = "se", N_v = "N"),
    transform = map(i = !!ind_i)
  ),
  input_d = target(map_merge_df(snp_df = snp_df, gwas_df = gwas_df),
                   transform = map(snp_df, gwas_df)),
  mmap_df = target(map_merge_map(input_d, mapd_df),
                   transform = map(input_d, mapd_df)),
  quh = target(map_local_quh_gen(mmap_df, "EUR"), transform = map(mmap_df)),
  est_df = target(map_local_rssp_est(quh), transform = map(quh)),
  full_est_df = target(rba(est_df), transform = combine(est_df))
)
