library(EigenH5)
library(ldmap)
library(dplyr)
data(ldetect_EUR)
inputf <- snakemake@input[["inputf"]]
cdf <- read_df_h5(inputf,"chrom_offset") %>% slice(1:22) 
cds <- as.integer(sum(cdf$datasize))

iv <- read_vector_h5(inputf, "snp/snp_struct", 1L:cds)
ldi <- snp_in_region(iv, ldetect_EUR)
rle2offset(ldi) %>% 
  rename(region_id=value) %>%
  mutate(offset=as.integer(offset)) %>% 
  saveRDS(snakemake@output[["offsetf"]])
