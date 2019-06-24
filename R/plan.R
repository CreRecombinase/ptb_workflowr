ld_df <- read_tsv(data_config$data$ldetect,col_types = c(chrom="i",start="i",stop="i",region_id="i"))

num_loci <- data_config$num_loci
db_df <-data_config$data$db
dcf <-data_config$data$anno
p_thresh <- data_config$p_thresh
all_feat <- str_replace(path_file(dcf),".bed.*","")

all_feat <- all_feat[str_detect(all_feat,"reproducible")|str_detect(all_feat,"seq",negate=T)]
plan <- drake_plan(
    gwas_df_ptb =  full_gwas_df(db_df,"beta","se","n",TRUE,nlines=data_config$nlines),
    p=calc_p(db_df),
    top_gwas_loc = group_by(gwas_df_ptb,region_id) %>% filter(abs(`z-stat`)==max(abs(`z-stat`))) %>% ungroup() %>%  arrange(desc(`z-stat`)),
    top_gwas_reg = top_gwas_loc %>% slice(1:num_loci),
    ngwas_i = mutate(gwas_df_ptb,SNP=1:n()),
    gr_df = make_range(gwas_df_ptb),
    ra=target(read_anno_r(feat_name,dcf=dcf),transform=map(feat_name=!!all_feat)),
    anno_r = target(anno_overlap_fun(input_range =ra,
                                     gr_df = gr_df,
                                     gw_df =ngwas_i),transform=map(ra)),
    s_feat = target(run_torus_Rdf(gw_df=ngwas_i,anno_df=anno_r),transform = map(anno_r)),
    all_feat_df = target(bind_rows(anno_r),transform=combine(anno_r)),
    allres_df=target(bind_results(s_feat),transform=combine(s_feat)),
    forward_feat_df = forward_reg_torus(res_df = allres_df,
                                        gw_df = ngwas_i,
                                        iternum = 5,
                                        p_thresh=p_thresh,
                                        anno_df = all_feat_df,
                                        prior = top_gwas_reg$region_id),
    susie_results = target(shim_susie(df = filter(forward_feat_df$prior,region_id==top_gwas_reg$region_id[x]),h_p = 0.25/p),transform=map(x=!!(1:num_loci)))
)
    # sub_split_top_gwas = semi_join(ngwas_i,top_gwas_reg) %>% split(.$region_id),
    # susie_i = map(as.character(top_gwas_reg$region_id),function(x,p,){
    #     h_p <- h/p
    #     idf <- inner_join(s_torus_pt$priors[[x]],sub_split_top_gwas[[x]],by=c("SNP","region_id"))
    #     sp <- h_p*nrow(idf)
    #     ret <- shim_susie(idf,R = Matrix::Diagonal(nrow(idf)),h_p = h_p)
    #     return(ret)
    # },p=nrow(gwas_df_ptb))


