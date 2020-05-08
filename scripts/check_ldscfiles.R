library(vroom)
  library(dplyr)
library(purrr)
library(readr)


  library(fs)

  ## save.image("wsl.RDS")
  ## stop()
  ## setwd("~/Dropbox/Repos/ldsc/workflow")
  ## load("wsl.RDS")

  ## yam
  ## l_file <- yaml::yaml.load_file("../workflow/annots.yaml")

  feat_list <- snakemake@params[["features"]]
  baseline_feat <- snakemake@params[["baseline_features"]]
  annot_name <- snakemake@params[["annot_name"]]
  l2dir <- snakemake@params[["L2"]]

  stopifnot(!is.null(feat_list),!is.null(baseline_feat),!is.null(annot_name),!is.null(l2dir))

file_df <- tibble::as_tibble(expand.grid(feature = feat_list, chrom = 1:22,stringsAsFactors = FALSE)) %>%
  dplyr::mutate(
           path = fs::path(l2dir,
                           paste0( feature, ".", chrom, ".l2.ldscore.gz")),
           annot_path = fs::path(l2dir,
                           paste0( feature, ".", chrom, ".annot.gz")),
           baseline_path = fs::path(l2dir,"baseline/", paste0("baselineLD.", chrom, ".l2.ldscore.gz")),
           baseline_annot_path = fs::path(l2dir,"baseline/", paste0("baselineLD.", chrom, ".annot.gz")),
           new_path = fs::path(l2dir,"new_baseline/", paste0(annot_name,".", chrom, ".l2.ldscore.gz")),
           new_annot_path = fs::path(l2dir,"new_baseline/", paste0(annot_name,".", chrom, "annot.gz"))
               )





stopifnot(all(fs::file_exists(c(file_df$path,file_df$baseline_path,file_df$annot_path,file_df$baseline_annot_path))))

spec <- cols(
  CHR = col_skip(),
  SNP = col_skip(),
  BP = col_skip(),
  L2 = col_double()
)
spec_base <- cols(
  CHR = col_double(),
  SNP = col_character(),
  BP = col_double(),
  baseL2 = col_double(),
  Coding_UCSCL2 = col_double(),
  Coding_UCSC.flanking.500L2 = col_double(),
  Conserved_LindbladTohL2 = col_double(),
  Conserved_LindbladToh.flanking.500L2 = col_double(),
  CTCF_HoffmanL2 = col_double(),
  CTCF_Hoffman.flanking.500L2 = col_double(),
  DGF_ENCODEL2 = col_double(),
  DGF_ENCODE.flanking.500L2 = col_double(),
  DHS_peaks_TrynkaL2 = col_double(),
  DHS_TrynkaL2 = col_double(),
  DHS_Trynka.flanking.500L2 = col_double(),
  Enhancer_AnderssonL2 = col_double(),
  Enhancer_Andersson.flanking.500L2 = col_double(),
  Enhancer_HoffmanL2 = col_double(),
  Enhancer_Hoffman.flanking.500L2 = col_double(),
  FetalDHS_TrynkaL2 = col_double(),
  FetalDHS_Trynka.flanking.500L2 = col_double(),
  H3K27ac_HniszL2 = col_double(),
  H3K27ac_Hnisz.flanking.500L2 = col_double(),
  H3K27ac_PGC2L2 = col_double(),
  H3K27ac_PGC2.flanking.500L2 = col_double(),
  H3K4me1_peaks_TrynkaL2 = col_double(),
  H3K4me1_TrynkaL2 = col_double(),
  H3K4me1_Trynka.flanking.500L2 = col_double(),
  H3K4me3_peaks_TrynkaL2 = col_double(),
  H3K4me3_TrynkaL2 = col_double(),
  H3K4me3_Trynka.flanking.500L2 = col_double(),
  H3K9ac_peaks_TrynkaL2 = col_double(),
  H3K9ac_TrynkaL2 = col_double(),
  H3K9ac_Trynka.flanking.500L2 = col_double(),
  Intron_UCSCL2 = col_double(),
  Intron_UCSC.flanking.500L2 = col_double(),
  PromoterFlanking_HoffmanL2 = col_double(),
  PromoterFlanking_Hoffman.flanking.500L2 = col_double(),
  Promoter_UCSCL2 = col_double(),
  Promoter_UCSC.flanking.500L2 = col_double(),
  Repressed_HoffmanL2 = col_double(),
  Repressed_Hoffman.flanking.500L2 = col_double(),
  SuperEnhancer_HniszL2 = col_double(),
  SuperEnhancer_Hnisz.flanking.500L2 = col_double(),
  TFBS_ENCODEL2 = col_double(),
  TFBS_ENCODE.flanking.500L2 = col_double(),
  Transcr_HoffmanL2 = col_double(),
  Transcr_Hoffman.flanking.500L2 = col_double(),
  TSS_HoffmanL2 = col_double(),
  TSS_Hoffman.flanking.500L2 = col_double(),
  UTR_3_UCSCL2 = col_double(),
  UTR_3_UCSC.flanking.500L2 = col_double(),
  UTR_5_UCSCL2 = col_double(),
  UTR_5_UCSC.flanking.500L2 = col_double(),
  WeakEnhancer_HoffmanL2 = col_double(),
  WeakEnhancer_Hoffman.flanking.500L2 = col_double(),
  GERP.NSL2 = col_double(),
  GERP.RSsup4L2 = col_double(),
  MAFbin1L2 = col_double(),
  MAFbin2L2 = col_double(),
  MAFbin3L2 = col_double(),
  MAFbin4L2 = col_double(),
  MAFbin5L2 = col_double(),
  MAFbin6L2 = col_double(),
  MAFbin7L2 = col_double(),
  MAFbin8L2 = col_double(),
  MAFbin9L2 = col_double(),
  MAFbin10L2 = col_double(),
  MAF_Adj_Predicted_Allele_AgeL2 = col_double(),
  MAF_Adj_LLD_AFRL2 = col_double(),
  Recomb_Rate_10kbL2 = col_double(),
  Nucleotide_Diversity_10kbL2 = col_double(),
  Backgrd_Selection_StatL2 = col_double(),
  CpG_Content_50kbL2 = col_double(),
  MAF_Adj_ASMCL2 = col_double(),
  GTEx_eQTL_MaxCPPL2 = col_double(),
  BLUEPRINT_H3K27acQTL_MaxCPPL2 = col_double(),
  BLUEPRINT_H3K4me1QTL_MaxCPPL2 = col_double(),
  BLUEPRINT_DNA_methylation_MaxCPPL2 = col_double(),
  synonymousL2 = col_double(),
  non_synonymousL2 = col_double(),
  Conserved_Vertebrate_phastCons46wayL2 = col_double(),
  Conserved_Vertebrate_phastCons46way.flanking.500L2 = col_double(),
  Conserved_Mammal_phastCons46wayL2 = col_double(),
  Conserved_Mammal_phastCons46way.flanking.500L2 = col_double(),
  Conserved_Primate_phastCons46wayL2 = col_double(),
  Conserved_Primate_phastCons46way.flanking.500L2 = col_double(),
  BivFlnkL2 = col_double(),
  BivFlnk.flanking.500L2 = col_double(),
  Human_Promoter_VillarL2 = col_double(),
  Human_Promoter_Villar.flanking.500L2 = col_double(),
  Human_Enhancer_VillarL2 = col_double(),
  Human_Enhancer_Villar.flanking.500L2 = col_double(),
  Ancient_Sequence_Age_Human_PromoterL2 = col_double(),
  Ancient_Sequence_Age_Human_Promoter.flanking.500L2 = col_double(),
  Ancient_Sequence_Age_Human_EnhancerL2 = col_double(),
  Ancient_Sequence_Age_Human_Enhancer.flanking.500L2 = col_double(),
  Human_Enhancer_Villar_Species_Enhancer_CountL2 = col_double(),
  Human_Promoter_Villar_ExACL2 = col_double(),
  Human_Promoter_Villar_ExAC.flanking.500L2 = col_double()
)



anno_cols <- cols(
  CHR = col_double(),
  BP = col_double(),
  SNP = col_character(),
  CM = col_double(),
  base = col_double(),
  Coding_UCSC = col_double(),
  Coding_UCSC.flanking.500 = col_double(),
  Conserved_LindbladToh = col_double(),
  Conserved_LindbladToh.flanking.500 = col_double(),
  CTCF_Hoffman = col_double(),
  CTCF_Hoffman.flanking.500 = col_double(),
  DGF_ENCODE = col_double(),
  DGF_ENCODE.flanking.500 = col_double(),
  DHS_peaks_Trynka = col_double(),
  DHS_Trynka = col_double(),
  DHS_Trynka.flanking.500 = col_double(),
  Enhancer_Andersson = col_double(),
  Enhancer_Andersson.flanking.500 = col_double(),
  Enhancer_Hoffman = col_double(),
  Enhancer_Hoffman.flanking.500 = col_double(),
  FetalDHS_Trynka = col_double(),
  FetalDHS_Trynka.flanking.500 = col_double(),
  H3K27ac_Hnisz = col_double(),
  H3K27ac_Hnisz.flanking.500 = col_double(),
  H3K27ac_PGC2 = col_double(),
  H3K27ac_PGC2.flanking.500 = col_double(),
  H3K4me1_peaks_Trynka = col_double(),
  H3K4me1_Trynka = col_double(),
  H3K4me1_Trynka.flanking.500 = col_double(),
  H3K4me3_peaks_Trynka = col_double(),
  H3K4me3_Trynka = col_double(),
  H3K4me3_Trynka.flanking.500 = col_double(),
  H3K9ac_peaks_Trynka = col_double(),
  H3K9ac_Trynka = col_double(),
  H3K9ac_Trynka.flanking.500 = col_double(),
  Intron_UCSC = col_double(),
  Intron_UCSC.flanking.500 = col_double(),
  PromoterFlanking_Hoffman = col_double(),
  PromoterFlanking_Hoffman.flanking.500 = col_double(),
  Promoter_UCSC = col_double(),
  Promoter_UCSC.flanking.500 = col_double(),
  Repressed_Hoffman = col_double(),
  Repressed_Hoffman.flanking.500 = col_double(),
  SuperEnhancer_Hnisz = col_double(),
  SuperEnhancer_Hnisz.flanking.500 = col_double(),
  TFBS_ENCODE = col_double(),
  TFBS_ENCODE.flanking.500 = col_double(),
  Transcr_Hoffman = col_double(),
  Transcr_Hoffman.flanking.500 = col_double(),
  TSS_Hoffman = col_double(),
  TSS_Hoffman.flanking.500 = col_double(),
  UTR_3_UCSC = col_double(),
  UTR_3_UCSC.flanking.500 = col_double(),
  UTR_5_UCSC = col_double(),
  UTR_5_UCSC.flanking.500 = col_double(),
  WeakEnhancer_Hoffman = col_double(),
  WeakEnhancer_Hoffman.flanking.500 = col_double(),
  GERP.NS = col_double(),
  GERP.RSsup4 = col_double(),
  MAFbin1 = col_double(),
  MAFbin2 = col_double(),
  MAFbin3 = col_double(),
  MAFbin4 = col_double(),
  MAFbin5 = col_double(),
  MAFbin6 = col_double(),
  MAFbin7 = col_double(),
  MAFbin8 = col_double(),
  MAFbin9 = col_double(),
  MAFbin10 = col_double(),
  MAF_Adj_Predicted_Allele_Age = col_double(),
  MAF_Adj_LLD_AFR = col_double(),
  Recomb_Rate_10kb = col_double(),
  Nucleotide_Diversity_10kb = col_double(),
  Backgrd_Selection_Stat = col_double(),
  CpG_Content_50kb = col_double(),
  MAF_Adj_ASMC = col_double(),
  GTEx_eQTL_MaxCPP = col_double(),
  BLUEPRINT_H3K27acQTL_MaxCPP = col_double(),
  BLUEPRINT_H3K4me1QTL_MaxCPP = col_double(),
  BLUEPRINT_DNA_methylation_MaxCPP = col_double(),
  synonymous = col_double(),
  non_synonymous = col_double(),
  Conserved_Vertebrate_phastCons46way = col_double(),
  Conserved_Vertebrate_phastCons46way.flanking.500 = col_double(),
  Conserved_Mammal_phastCons46way = col_double(),
  Conserved_Mammal_phastCons46way.flanking.500 = col_double(),
  Conserved_Primate_phastCons46way = col_double(),
  Conserved_Primate_phastCons46way.flanking.500 = col_double(),
  BivFlnk = col_double(),
  BivFlnk.flanking.500 = col_double(),
  Human_Promoter_Villar = col_double(),
  Human_Promoter_Villar.flanking.500 = col_double(),
  Human_Enhancer_Villar = col_double(),
  Human_Enhancer_Villar.flanking.500 = col_double(),
  Ancient_Sequence_Age_Human_Promoter = col_double(),
  Ancient_Sequence_Age_Human_Promoter.flanking.500 = col_double(),
  Ancient_Sequence_Age_Human_Enhancer = col_double(),
  Ancient_Sequence_Age_Human_Enhancer.flanking.500 = col_double(),
  Human_Enhancer_Villar_Species_Enhancer_Count = col_double(),
  Human_Promoter_Villar_ExAC = col_double(),
  Human_Promoter_Villar_ExAC.flanking.500 = col_double()
)

file_df <- nest(file_df,feat_data = c(feature,path,annot_path))
  modify_cols <- function(cols,old,new) {
    cols$cols <- set_names(cols$cols, ~dplyr::if_else(.== old,new, .))
    return(cols)
  }


  feat_fun <- function(f,feat) {
    return(tibble::tibble(!!feat := scan(f,what = numeric(),skip = 1)))
    }


  pwalk(file_df,function(baseline_path, chrom, feat_data, new_path, new_annot_path,baseline_annot_path) {
    cat("Now on",chrom,"\n")
    keep_cols <- c("CHR","SNP","BP","baseL2",baseline_feat)

    bldcols <- names(spec_base$cols)
    bacols <- names(anno_cols$cols)
    stopifnot(all(keep_cols %in%  names(spec_base$cols)))
    stopifnot(all(keep_cols %in%  names(anno_cols$cols)))
    bad_good_cols <- keep_cols[!keep_cols %in%  names(anno_cols$cols)]
  
    bad_cols <- names(spec_base$cols)[!(names(spec_base$cols) %in% keep_cols)]
    bad_anno_cols <- str_replace(bad_cols,"L2$","")

    new_base <- spec_base
    new_anno_spec <- anno_cols
    for (bc in seq_along(bad_cols)) {
      new_base$cols[[bad_cols[bc]]] <- col_skip()
      new_anno_spec$cols[[bad_anno_cols[bc]]] <- col_skip()
    }
  
    feat_data <- mutate(feat_data,col_spec = map(feature, ~modify_cols(spec, "L2" , .x)))
    annot_df <- map2_dfc(feat_data$annot_path,feat_data$feature, ~feat_fun(.x, .y))
    o_anno_path <- bind_cols(vroom::vroom(baseline_annot_path,col_types = new_anno_spec,delim = "\t"),annot_df)
    tannot <- read_delim(feat_data$annot_path[1],delim = "\t")
    lddf <- pmap_dfc(feat_data,function(feature,path,col_spec, ...) {
      vroom::vroom(path,col_names = names(col_spec$cols),col_types = col_spec,delim = "\t",skip = 1L)})
    baseline_df <- bind_cols(vroom::vroom(baseline_path,delim = "\t",col_types = new_base),lddf)
    vroom::vroom_write(baseline_df,new_path,delim = "\t")
  })




om(baseline_path,delim = "\t",col_types = new_base),lddf)
    vroom::vroom_write(baseline_df,new_path,delim = "\t")
  })
