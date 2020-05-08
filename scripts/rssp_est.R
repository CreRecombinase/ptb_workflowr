library(dplyr)
library(ldmap)
library(ldshrink)
library(EigenH5)
library(purrr)

input_f <- snakemake@input[["h5f"]]
sumstat_f <- snakemake@input[["gwasf"]]
sumstat_df <- vroom::vroom(sumstat_f) %>%
  mutate(SNP = rsid2int(SNP))

read_rsid_h5 <- function(file, ldmr_id){
  tibble(ldmr=as.integer(ldmr_id),rsid=read_vector_h5(file, paste0(ldmr_id, "/rsid"))) %>%
    mutate(ld_id=1:n())
}

all_reg <- map2_dfr(input_f, 1:22, function(x, y){
  map_dfr(ls_h5(x),~read_rsid_h5(x,.x)) %>% mutate(chrom=y)
})

squh_h5 <- function(ldmr,Z,ld_f,rsid){
    Q <- read_matrix_h5(ld_f,fs::path("Q",as.character(ldmr[1])))
    tZ <- Z[as.character(rsid),,drop=FALSE]
    return(RSSp::convert_quh(tZ,Q))
}


ld_df <- inner_join(all_reg,sumstat_df,by=c("rsid"="SNP"))
ld_dfl <- split(ld_df,ld_df$ldmr)

read_r <- function(tdf){
  R <- read_matrix_h5(input_f[unique(tdf$chrom)],paste0(tdf$ldmr[1],"/R"))[tdf$ld_id,tdf$ld_id,drop=FALSE]
  ldvr <- eigen(R)
  quh <- RSSp::convert_quh(tdf$Z,ldvr$vectors)
  tibble(quh=quh,D=ldvr$values)
}

rssp_df <- map_df(ld_dfl,read_r)

res_df <- RSSp::RSSp_estimate(rssp_df$quh,rssp_df$D,mean(sumstat_df$N))
saveRDS(res_df,snakemake@output[["est_rdsf"]])
