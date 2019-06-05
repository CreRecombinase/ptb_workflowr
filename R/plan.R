ld_df <- read_tsv(ldetect_file,col_types = c(chrom="i",start="i",stop="i",region_id="i"))
t_ld_df <- mutate(ld_df) %>% select(-region_id)
ld_df_l <- map(parallel::splitIndices(nrow(ld_df),50),~slice(ld_df,.x))
ind_i <- seq_along(ld_df_l)


all_feat <- fs::path_ext_remove(fs::path_file(fs::dir_ls(data_config$anno_dir,glob="*bed")))
snp_feat <- fs::path_ext_remove(fs::path_file(fs::dir_ls(data_config$anno_dir,glob="*zscore.gz")))

b_eqtl <- all_feat[!str_detect(all_feat,"eQTL")]
plan <- drake_plan(
    gwas_df_ptb =  full_gwas_df("beta","se","N"),
    gr_df = target(make_range(x),transform = map(x= c(gwas_df_ptb))),
    readr::write_tsv(dplyr::select(gwas_df_ptb,SNP,locus = region_id,`z-stat`),file_out("gwas_ptb_file.tsv.gz")),
    w_eqtl = target(
        do_torus(c("eQTL_0.05_FDR",f),
                 gr_df=gr_df_gwas_df_ptb,
                 gw_df = gwas_df_ptb,
                 gwas_file = file_in("gwas_ptb_file.tsv.gz")),
        transform = map(f=!!b_eqtl)),
    cust_t = do_torus_p(
        c("E8_TCM_D_48h",
          "eQTL_0.05_FDR",
          "FOSL2",
          "h3k27ac-final-C-48h",
          "FOXO1_DeMayo",
          "Repressed_Hoffman",
          "POLII_DeMayo"),
                 gr_df=gr_df_gwas_df_ptb,
                 gw_df = gwas_df_ptb,
                 gwas_file = file_in("gwas_ptb_file.tsv.gz")),
    all_we = target(bind_rows(w_eqtl),transform=combine(w_eqtl)),
    all_wel = target(list(w_eqtl),transform=combine(w_eqtl)),
    all_wel_df = imap_dfr(all_wel,~mutate(.x,id=.y))
)



#
gf_n <- tribble(
    ~term, ~name,
    "Glut_FDR_0.05.pgc2.zscore.1",                      "iN_Glut_ASoC",
    "NPC_FDR_0.05.pgc2.zscore.1",                       "NPC_ASoC",
    "GA_all_peaks.narrowPeak.cleaned.hg19.merged.1",  "iN_GABA_OCR",
    "DN_all_peaks.narrowPeak.cleaned.hg19.merged.1",  "iN_DA_OCR",
    "CN_all_peaks.narrowPeak.cleaned.hg19.merged.1",  "iN_Glut_OCR",
    "ips_all_peaks.narrowPeak.cleaned.hg19.merged.1",  "iPSC_OCR",
    "NSC_all_peaks.narrowPeak.cleaned.hg19.merged.1",  "NPC_OCR",
     "Conserved_LindbladToh.1"                      ,  "Conserved",
     "Coding_UCSC.1"                      ,  "Coding",
    "Promoter_UCSC.1"                      ,  "Promoter",
)



# f_all_res_scz <- bind_rows(all_res_scz,snp_all_res_scz) %>% filter(term!="Intercept") %>%
#     inner_join(gf_n) %>%
#     mutate(estimate=estimate/log(2),low=low/log(2),high=high/log(2))
#
# # iPSC_OCR: ips
# # NPC_OCR: NSC
# # iN_Glut_OCR: CN
# # iN_DA_OCR: DN
# # iN_GABA_OCR: GA
#     # gr_df = target(dplyr::select(x, chrom, start = pos) %>%
#     # mutate(chrom = paste0("chr", chrom),
#     #        start = as.integer(start),
#     #        end = start + 1L) %>%
#     # GenomicRanges::makeGRangesFromDataFrame(),
#     # gwas_r =  gwas_range(),
#
#     # anno_hic = read_anno_r("hic_harmonizome_interacting_DT1_dTL4_D_48h"),
#     # hic_gw_df = anno_overlap_fun(anno_hic,gr_df,gwas_df,"hic_harmonizome_interacting_DT1_dTL4_D_48h"),
#
#     # eqtl_res = run_torus(file_in("gwas_file.tsv.gz"),file_in("eqtl_gwas_file.tsv.gz")),
#     # ut_res = run_torus(file_in("gwas_file.tsv.gz"),file_in("ut_eqtl.tsv.gz")),
#     # hic_res = run_torus(file_in("gwas_file.tsv.gz"),file_in("hic_t.tsv.gz")),
#     # big_hic_res = run_torus(file_in("gwas_file.tsv.gz"),file_in("big_hic_t.tsv.gz")),
#     # fat_res = run_torus(file_in("gwas_file.tsv.gz"),file_in("fat_eqtl.tsv.gz")),
# # big_eqtl_df = vroom::vroom("~/Downloads/GTEx_Analysis_v7_meta_tmp.tsv",col_types = cols_only(
# #   RSID = col_character(),
# #   `#STUDY` = col_double(),
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
