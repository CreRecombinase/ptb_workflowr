num_loci <- data_config$num_loci
max_size <- data_config$max_snp
min_snp <- data_config$min_snp
db_df <-data_config$data$db
dcf <-data_config$data$anno
p_thresh <- data_config$p_thresh
#all_feat <- str_replace(path_file(dcf),".bed.*","")
h <- 0.25
geno_f <- data_config$data$ldp
input_feat_df <- tibble::tibble(term=all_feat,estimate=NA_real_,low=NA_real_,high=NA_real_,sd=NA_real_,z=NA_real_,p=0,lik=-Inf)
saf <- seq_along(all_feat)
best_terms <- c("atac-seq-pooled-DSC1-dec-ATAC",
                "chip-seq-dec_up-H3K27ac",
                "chip-seq-dec_up-H3K4me1",
                "chip-seq-pooled-DSC2-ctr-H3K4me1",
                "eQTL_0.05_FDR")
all_feat <- unique(c(best_terms,all_feat))

#Intercept    -11.330       -11.352    -11.307
#atac-seq-pooled-DSC1-dec-ATAC.1     -4.089        -6.320     -1.859
#chip-seq-dec_up-H3K27ac.1      2.784         1.453      4.116
#chip-seq-dec_up-H3K4me1.1     -2.539        -3.677     -1.400
#chip-seq-pooled-DSC1-ctr-H3K27ac.1     -0.883        -1.499     -0.267
#chip-seq-pooled-DSC2-ctr-H3K4me1.1      3.143         2.253      4.032
#eQTL_0.05_FDR.1      6.460         1.590     11.329
p <- 14991823
plan <- drake_plan(
    sgwas_df_ptb =  target(full_gwas_df(db_df = db_df,
                                        beta_v = "beta",
                                        se_v = "se",
                                        N_v = "n",
                                        p_v = "pval",
                                        keep_bh_se = TRUE,
                                        keep_allele = TRUE,
                                        nlines=data_config$nlines
                                        ),trigger=trigger(change=c(db_df,data_config$nlines),command = T,file = F,depend = T)),
    pre_gwas_df_ptb = assign_reg_df(sgwas_df_ptb,data_config$data$ld_df,max_snp=max_size,min_snp=min_snp),
    gwas_df_ptb = merge_snp_f(file_in(data_config$data$ldp),gwas_df = pre_gwas_df_ptb),
    top_gwas_loc = group_by(gwas_df_ptb,region_id) %>%
        filter(abs(`z-stat`)==max(abs(`z-stat`))) %>%
        ungroup() %>%
        arrange(desc(`z-stat`)) %>% mutate(pz=percent_rank(`z-stat`)),
    top_gwas_reg = top_gwas_loc %>%
        slice(1:num_loci),
    slice_gw_df = filter(top_gwas_loc,pz>0.5) %>% dplyr::select(region_id) %>% inner_join(gwas_df_ptb),
    gr_df = make_range(slice_gw_df),
    ra=target(read_anno_r(feat_name,dcf=dcf),transform=map(feat_name=!!best_terms)),
    anno_r = target(anno_overlap_fun(input_range =ra,
                                     gr_df = gr_df,
                                     gw_df =slice_gw_df),transform=map(ra)),
    taf =  target(write_anno(anno_r,p),transform = map(anno_r)),
    naf = write_anno(p = p),
    full_anno_df = target(bind_rows(anno_r),transform=combine(anno_r)),
    faf = write_anno(dplyr::filter(full_anno_df,feature %in% best_terms)),
    gf = write_gwas(slice_gw_df),
    ind_results =  target(run_torus_cmd(gf = gf,af = taf),transform = map(taf)),
    null_results = run_torus_cmd(gf = gf,af = naf,torus_p = top_gwas_reg$region_id),
    mult_results = run_torus_cmd(gf = gf,af = faf,torus_p = top_gwas_reg$region_id),
    split_gw_df  = semi_join(slice_gw_df,dplyr::select(top_gwas_reg,region_id)) %>% split(.$region_id),
    prior_r = purrr::map(top_gwas_reg$region_id,~dplyr::inner_join(mult_results$priors[[.x]],split_gw_df[[.x]])),
    null_r =  purrr::map(top_gwas_reg$region_id,~dplyr::inner_join(null_results$priors[[.x]],split_gw_df[[.x]])),
    susie_res = target(shim_susie(prior_r[[ix]],10,h_p=h/p,geno_f=geno_f),transform=map(ix = !!seq_len(data_config$num_loci))),
    nullsusie_res = target(shim_susie(null_r[[ix]],10,h_p=h/p,geno_f=geno_f),transform=map(ix = !!seq_len(data_config$num_loci))),
    f_susie = target(merge_susie(prior_r[[ix]],susie_res),transform=map(ix=!!seq_len(data_config$num_loci),susie_res)),
    fnull_susie = target(merge_susie(null_r[[ix]],nullsusie_res),transform=map(ix=!!seq_len(data_config$num_loci),nullsusie_res)),



)
    # sub_split_top_gwas = semi_join(gwas_df_ptb,top_gwas_reg) %>% split(.$region_id),
    # susie_i = map(as.character(top_gwas_reg$region_id),function(x,p,){
    #     h_p <- h/p
    #     idf <- inner_join(s_torus_pt$priors[[x]],sub_split_top_gwas[[x]],by=c("SNP","region_id"))
    #     sp <- h_p*nrow(idf)
    #     ret <- shim_susie(idf,R = Matrix::Diagonal(nrow(idf)),h_p = h_p)
    #     return(ret)
    # },p=nrow(gwas_df_ptb))


