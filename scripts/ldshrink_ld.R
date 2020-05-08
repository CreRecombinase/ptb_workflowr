library(dplyr)
library(ldmap)
library(ldshrink)
library(EigenH5)
shrink <- snakemake@params[["shrink"]]
if(is.null(shrink)){
  doshrink <- TRUE
}else{
  doshrink <- shrink=="shrink"
}
bim_df <- read_plink_bim(snakemake@input[["bimf"]]) %>% 
  mutate(snp_id = 1:n(),
         ldmr = snp_overlap_region(snp_struct, ldetect_EUR),
         rsid=rsid2int(rsid))
fam_df <- read_plink_fam(snakemake@input[["famf"]])
N <- nrow(fam_df)
bim_l <- split(bim_df, bim_df$ldmr)
purrr::walk(bim_l, function(df){
  gl <- read_plink_bed(snakemake@input[["bedf"]], subset = df$snp_id, N = N)
  Xm <- gt2matrix(gl)
  if(!doshrink){
    R <- stats::cor(Xm, use = "complete.obs")
  }else{
    R <- ldshrink::ldshrink(Xm, df$map, isGeno = TRUE)
  }
  ldmr_id <- as.character(unique(df$ldmr))
  write_matrix_h5(R, snakemake@output[["h5f"]], paste0(ldmr_id, "/R"))
  write_vector_h5(df$snp_id, snakemake@output[["h5f"]], paste0(ldmr_id, "/snp_id"))
  write_vector_h5(df$rsid, snakemake@output[["h5f"]], paste0(ldmr_id, "/rsid"))
})
