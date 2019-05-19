ld_df <- read_tsv(ldetect_file)
t_ld_df <- mutate(ld_df) %>% select(-region_id)
ld_df_l <- map(parallel::splitIndices(nrow(ld_df),50),~slice(ld_df,.x))
ind_i <- seq_along(ld_df_l)
#http://shiny.stephenslab.uchicago.edu/gaow/finemapping_summary_stats_example.tar.gz


all_feat <- fs::path_ext_remove(fs::path_file(fs::dir_ls(data_config$anno_dir,glob="*bed")))

plan <- drake_plan(
    gwas_df_ptb =  full_gwas_df("beta","se","N"),
    gwas_df_scz = read_scz(),
    gr_df = target(make_range(x),transform = map(x= c(gwas_df_scz,gwas_df_ptb))),
    readr::write_tsv(dplyr::select(gwas_df_ptb,SNP,locus = region_id,`z-stat`),file_out("gwas_ptb_file.tsv.gz")),
    readr::write_tsv(dplyr::select(gwas_df_scz,SNP,locus = region_id,`z-stat`),file_out("gwas_scz_file.tsv.gz")),
    res_df_ptb = target(
        do_torus(f,
                 gr_df=gr_df_gwas_df_ptb,
                 gw_df = gwas_df_ptb,
                 gwas_file = file_in("gwas_ptb_file.tsv.gz")),
        transform = map(f=!!all_feat)),
    res_df_scz = target(
        do_torus(f,
                 gr_df=gr_df_gwas_df_scz,
                 gw_df = gwas_df_scz,
                 gwas_file = file_in("gwas_scz_file.tsv.gz")),
        transform = map(f=!!all_feat)),
    all_res_scz = target(bind_rows(res_df_scz),transform=combine(res_df_scz)),
    all_res_ptb = target(bind_rows(res_df_ptb),transform=combine(res_df_ptb))



)

    # gr_df = target(dplyr::select(x, chrom, start = pos) %>%
    # mutate(chrom = paste0("chr", chrom),
    #        start = as.integer(start),
    #        end = start + 1L) %>%
    # GenomicRanges::makeGRangesFromDataFrame(),
    # gwas_r =  gwas_range(),

    # anno_hic = read_anno_r("hic_harmonizome_interacting_DT1_dTL4_D_48h"),
    # hic_gw_df = anno_overlap_fun(anno_hic,gr_df,gwas_df,"hic_harmonizome_interacting_DT1_dTL4_D_48h"),

    # eqtl_res = run_torus(file_in("gwas_file.tsv.gz"),file_in("eqtl_gwas_file.tsv.gz")),
    # ut_res = run_torus(file_in("gwas_file.tsv.gz"),file_in("ut_eqtl.tsv.gz")),
    # hic_res = run_torus(file_in("gwas_file.tsv.gz"),file_in("hic_t.tsv.gz")),
    # big_hic_res = run_torus(file_in("gwas_file.tsv.gz"),file_in("big_hic_t.tsv.gz")),
    # fat_res = run_torus(file_in("gwas_file.tsv.gz"),file_in("fat_eqtl.tsv.gz")),
# big_eqtl_df = vroom::vroom("~/Downloads/GTEx_Analysis_v7_meta_tmp.tsv",col_types = cols_only(
#   RSID = col_character(),
#   `#STUDY` = col_double(),
#   PVALUE_FE = col_double(),
#   BETA_FE = col_double(),
#   STD_FE = col_double(),
#   PVALUE_RE = col_double(),
#   BETA_RE = col_double(),
#   STD_RE = col_double(),
#   PVALUE_RE2 = col_double(),
#   STAT1_RE2 = col_double(),
#   STAT2_RE2 = col_double(),
#   I_SQUARE = col_double(),
#   Q = col_double(),
#   PVALUE_Q = col_double(),
#   TAU_SQUARE = col_double()
# ),delim="\t") %>% filter(`#STUDY`>1),
## plan <- drake_plan(
##   snp_df = target(
##     map_df_reg("EUR", ld_df_l[[i]]),
##     transform = map(i = !!ind_i)),
##   mapd_df = target(map_read_map(ld_df_l[[i]]),
##                    transform = map(i = !!ind_i)),
##   gwas_df = target(
##     map_snp_reg(ld_df_l[[i]], beta_v = "beta", se_v = "se", N_v = "N"),
##     transform = map(i = !!ind_i)
##   ),
##   input_d = target(map_merge_df(snp_df = snp_df, gwas_df = gwas_df),
##                    transform = map(snp_df, gwas_df)),
##   mmap_df = target(map_merge_map(input_d, mapd_df),
##                    transform = map(input_d, mapd_df)),
##   quh = target(map_local_quh_gen(mmap_df, "EUR"), transform = map(mmap_df)),
##   est_df = target(map_local_rssp_est(quh), transform = map(quh)),
##   full_est_df = target(rba(est_df), transform = combine(est_df))
## )
