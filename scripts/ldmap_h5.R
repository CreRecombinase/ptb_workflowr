library(fs)
library(EigenH5)
library(ldmap)
library(dplyr)
library(purrr)

iff <- snakemake@input[["bimlist"]]
off <- snakemake@output[["h5"]]

bim_df <- read_plink_bim(snakemake@input[["bimlist"]])
fam_df <- read_plink_fam(snakemake@input[["famlist"]])
N <- as.integer(nrow(fam_df))
p <- as.integer(nrow(bim_df))
EigenH5::create_matrix_h5(off,"dosage",data=numeric(),dim=c(N,p))
bim_df <- bim_df %>% dplyr::mutate(index = 1:dplyr::n())
bim_df %>% EigenH5::write_df_h5(off,"snp")
i_l <- split(bim_df$index,ggplot2::cut_interval(x=seq_len(p),length=50000))
walk(i_l,function(subset){
  gl <- gt2matrix(read_plink_bed(snakemake@input[["bedlist"]],subset=subset,N=N))
  EigenH5::write_matrix_h5(gl, off, "dosage", subset_cols=subset)
})
