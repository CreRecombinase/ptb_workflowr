library(dplyr)
vroom::vroom(snakemake@input[["vecf"]]) %>% 
  dplyr::transmute(SNP=ID,A1,A2=REF,N=OBS_CT,Z=BETA/SE) %>% 
vroom::vroom_write(snakemake@output[["tempf"]],delim="\t")
