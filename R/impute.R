
snp_f <- fs::path(data_config$data$scratch,"SNPloc_all.txt")
snp_bedf <- fs::path(data_config$data$scratch,"SNPloc_all.bed")

fa_f <- fs::path(data_config$data$scratch,"hg19.fa.gz")
fa_i <- Rsamtools::indexFa(fa_f)
idx <- Rsamtools::scanFaIndex(fa_i)
snp_df <- readr::read_tsv(snp_f)
mutate(snp_df,p2=pos+1) %>% dplyr::select(chrm_snp,pos,p2) %>% write_tsv(path = snp_bedf,col_names = FALSE)
