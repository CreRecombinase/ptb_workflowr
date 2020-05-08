library(dplyr)
vroom::vroom(snakemake@input[["freqf"]]) %>% filter(ALT_FREQS>as.numeric(snakemake@params[["snp_freq"]])) %>% 
select(ID) %>% 
readr::write_tsv(snakemake@output[["snp_list"]],col_names=FALSE)
