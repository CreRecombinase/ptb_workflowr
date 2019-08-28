num_loci <- data_config$num_loci
max_size <- data_config$max_snp
min_snp <- data_config$min_snp
db_df <- data_config$data$db
dcf <- data_config$data$anno
p_thresh <- data_config$p_thresh

h <- 0.25
geno_f <- data_config$data$ldp

p <- 14991823
plan <- drake_plan(
    sgwas_df_ptb =  target(
        full_gwas_df(db_df = db_df,
                     beta_v = "beta",
                     se_v = "se",
                     N_v = "n",
                     p_v = "pval",
                     keep_bh_se = TRUE,
                     keep_allele = TRUE,
                     nlines = data_config$nlines
                     ),
        trigger = trigger(
            change = c(db_df, data_config$nlines),
            command = T,
            file = F,
            depend = T)),
    pre_gwas_df_ptb = assign_reg_df(sgwas_df_ptb,
                                    data_config$data$ld_df,
                                    max_snp = max_size,
                                    min_snp = min_snp),
    gwas_df_ptb = merge_snp_f(
        file_in(data_config$data$ldp),
        gwas_df = pre_gwas_df_ptb
    ),
    top_gwas_loc = group_by(gwas_df_ptb, region_id) %>%
        filter(p == min(p)) %>%
        ungroup() %>%
        arrange(p) %>% mutate(pz = percent_rank(p)),
    top_gwas_reg = top_gwas_loc %>%
        slice(1:num_loci),
    gr_df = make_range(gwas_df_ptb),
    ra = target(read_anno_r(feat_name, dcf = dcf),
                transform = map(feat_name = !!all_feat)),
    anno_r = target(anno_overlap_fun(input_range = ra,
                                     gr_df = gr_df,
                                     gw_df = gwas_df_ptb),
                    transform = map(ra)),
    target(write_anno(anno_r, af = file_out(pf)),
           transform = map(anno_r, pf = !!af)),
    full_anno_df = target(bind_rows(anno_r),
                          transform = combine(anno_r)),
    target(write_anno(
        dplyr::filter(full_anno_df,
                      feature %in% bt), af = file_out(f)),
        transform = map(bt = !!model_df$features, f = !!(model_df$file))),
    write_gwas(gwas_df_ptb, gf = file_out(gf)),
    ind_results =  target(
        run_torus_cmd(gf = file_in(gf),
                      af = file_in(taf)),
        transform = map(taf = !!af)),
    mix_results = target(
        run_torus_cmd(gf = file_in(gf),
                      af =  file_in(f),
                      torus_p = top_gwas_reg$region_id),
        transform = map(f = !!model_df$file)),
    split_gw_df  = semi_join(gwas_df_ptb,
                             dplyr::select(top_gwas_reg,
                                           region_id)) %>% split(.$region_id),
    prior_r = target(
        purrr::map(top_gwas_reg$region_id,
                   ~dplyr::inner_join(
                               mix_results$priors[[.x]],
                               split_gw_df[[.x]])),
        transform = map(mix_results)),
    susie_res = target(
        shim_susie(
            prior_r[[ix]],
            L = tL,
            h_p = h / p,
            geno_f = geno_f),
        transform = cross(
            prior_r,
            ix = !!seq_len(data_config$num_loci),
            tL = c(1L,2L,3L))),
)
