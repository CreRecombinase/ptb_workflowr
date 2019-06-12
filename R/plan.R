ld_df <- read_tsv(ldetect_file,col_types = c(chrom="i",start="i",stop="i",region_id="i"))
t_ld_df <- mutate(ld_df) %>% select(-region_id)
ld_df_l <- map(parallel::splitIndices(nrow(ld_df),50),~slice(ld_df,.x))
ind_i <- seq_along(ld_df_l)

data_config$anno_dir <- "/home/nwknoblauch/Dropbox/scratch/ptb_scratch/new_bed/"
all_feat <-fs::path_ext_remove(fs::path_ext_remove(fs::path_file(fs::dir_ls(data_config$anno_dir,glob="*bed*",type = "file"))))
all_feat <- c(all_feat,"eQTL_0.05_FDR","Repressed_Hoffman")

b_eqtl <- all_feat[!str_detect(all_feat,"eQTL")]
plan <- drake_plan(
    gwas_df_ptb =  full_gwas_df("beta","se","N",TRUE),
    top_gwas_reg = group_by(gwas_df_ptb,region_id) %>% summarise(top_z=max(abs(`z-stat`))) %>% arrange(desc(top_z)) %>% slice(1:50),
    ngwas_i = mutate(gwas_df_ptb,SNP=1:n()),
    gr_df = make_range(gwas_df_ptb),
    s_feat = target(anno_overlap_fun(input_range = list(a=read_anno_r(feat_name)),
                   gr_df = gr_df,
                   gw_df =ngwas_i,
                   name = feat_name),transform = map(feat_name=!!all_feat)),
    all_feat_df = target(bind_annotations(s_feat),transform=combine(s_feat)),
    feat_res =target(torusR(ngwas_i,list(s_feat)),transform=map(s_feat)),
    allres_df=target(bind_rows(feat_res),transform=combine(feat_res)),
    forward_feat_df = forward_reg_torus(res_df = allres_df,gw_df = ngwas_i,anno_df = all_feat_df),
    h_p=0.245/nrow(ngwas_i),
    susie_results = target(shim_susie(forward_feat_df$prior_l[[x]],h_p = h_p),transform=map(x=!!(1:50)))
)
    # sub_split_top_gwas = semi_join(ngwas_i,top_gwas_reg) %>% split(.$region_id),
    # susie_i = map(as.character(top_gwas_reg$region_id),function(x,p,h=.245){
    #     h_p <- h/p
    #     idf <- inner_join(s_torus_pt$priors[[x]],sub_split_top_gwas[[x]],by=c("SNP","region_id"))
    #     sp <- h_p*nrow(idf)
    #     ret <- shim_susie(idf,R = Matrix::Diagonal(nrow(idf)),h_p = h_p)
    #     return(ret)
    # },p=nrow(gwas_df_ptb))


