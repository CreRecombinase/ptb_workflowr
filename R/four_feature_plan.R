ld_df <- read_tsv(data_config$data$ldetect,col_types = c(chrom="i",start="i",stop="i",region_id="i"))

num_loci <- data_config$num_loci
max_size <- data_config$max_snp
min_snp <- data_config$min_snp
db_df <-data_config$data$db
dcf <-data_config$data$anno
p_thresh <- data_config$p_thresh

plan <- drake_plan(
  sgwas_df_ptb =  target(full_gwas_df(db_df,"beta","se","n",TRUE,nlines=data_config$nlines)),
  gwas_df_ptb = assign_reg_df(sgwas_df_ptb,ld_df,max_snp=max_size,min_snp=min_snp),
  p=calc_p(db_df),
  top_gwas_loc = group_by(gwas_df_ptb,region_id) %>% filter(abs(`z-stat`)==max(abs(`z-stat`))) %>% ungroup() %>%  arrange(desc(`z-stat`)),
  top_gwas_reg = top_gwas_loc %>% slice(1:num_loci),
  gr_df = make_range(gwas_df_ptb),
  ra=target(read_anno_r(feat_name,dcf=dcf),transform=map(feat_name=!!all_feat)),
  anno_r = target(anno_overlap_fun(input_range =ra,
                                   gr_df = gr_df,
                                   gw_df =gwas_df_ptb),transform=map(ra)),
  s_feat = target(run_torus_Rdf(gw_df=gwas_df_ptb,anno_df=anno_r),transform = map(anno_r)),
  all_feat_df = target(bind_rows(anno_r),transform=combine(anno_r)),
    feat4 = forward_reduce(f_feat = feat3,term_list = all_feat,lik_vec = alik4),
    prior_r = pr_torus(gw_df = gwas_df_ptb,anno_df = all_feat_df,feat_v = feat4,prior = top_gwas_reg$region_id),
    susie_res = target(shim_susie(df = prior_r$prior[[tix]],h_p=0.25/p),transform=map(tix=!!(1:num_loci)))
)
