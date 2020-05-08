library(ldmap)
library(dplyr)
fam_df <- read_plink_fam(snakemake@input[["famf"]])
ind_df <- read.table(snakemake@input[["grm_id"]],header=FALSE,stringsAsFactors=FALSE) %>% 
rename(fid=V1,iid=V2) %>% 
mutate(fid=as.character(fid),iid=as.character(iid))
rest_df <- anti_join(fam_df,ind_df) %>% sample_n(10000,replace=F)
write_plink_fam(snakemake@output[["sub_f"]])
