

ldetect_file <- fs::path_expand(fs::path(data_config$LDETECT_DIR,"fourier_ls-all.tsv.gz"))
#kg_dir <-fs::path_expand(fs::path(data_config$KG_DIRECTORY,"h5/EUR/"))
input_db_f <- fs::path_expand(fs::path(data_config$KG_DIRECTORY,"ga_gwas.db"))



stopifnot(all(file.exists(c(
                          ldetect_file,
                          input_db_f))))
