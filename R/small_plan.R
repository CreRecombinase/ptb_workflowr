ld_df <- read_tsv(ldetect_file,col_types = c(chrom="i",start="i",stop="i",region_id="i"))
ld_df_l <- map(parallel::splitIndices(nrow(ld_df),50),~slice(ld_df,.x))
ind_i <- seq_along(ld_df_l)


#all_feat <- fs::path_ext_remove(fs::path_file(fs::dir_ls(data_config$anno_dir,glob="*bed")))
data_config$anno_dir <- "/home/nwknoblauch/Dropbox/scratch/ptb_scratch/new_bed/"
all_feat <-fs::path_ext_remove(fs::path_ext_remove(fs::path_file(fs::dir_ls(data_config$anno_dir,glob="*bed*",type = "file"))))
all_feat <- c(all_feat,"eQTL_0.05_FDR","Repressed_Hoffman","hic_all_interacting_DT1_dTL4_D_48h")
all_feat <- all_feat[str_detect(all_feat,"reproducible")|str_detect(all_feat,"seq",negate=T)]


plan <- drake_plan(
    little_gwas = read_ptb_db(input_db_f) %>%
      head(5000) %>%
      collect() %>% mutate(`z-stat`=beta/se) %>%
      assign_reg_df(ld_df),
    little_gr_df = make_range(little_gwas),
    little_ngwas_i = mutate(little_gwas,SNP=1:n()),
     readr::write_tsv( dplyr::select(little_ngwas_i,SNP,locus = region_id,`z-stat`),file_write(file_out("little_gwas_i.tsv.zstd"),"zstd")),
    s_torus_pt = do_torus_p(
        feat_name = "Conserved_LindbladToh",
        gr_df=little_gr_df,
        gw_df = little_ngwas_i,
        gwas_file = file_in("little_gwas_i.tsv.zstd")),
)

