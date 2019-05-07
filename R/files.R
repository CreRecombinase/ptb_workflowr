

map_file = fs::path_expand(fs::path(data_config$KG_MAPDIR,"interpolated_EUR.h5"))
ldetect_file <- fs::path_expand(fs::path(data_config$LDETECT_DIR,"fourier_ls-all.tsv.gz"))
kg_dir <-fs::path_expand(fs::path(data_config$KG_DIRECTORY,"h5/EUR/"))
input_db_f <- fs::path_expand(fs::path(data_config$KG_DIRECTORY,"ga_gwas.db"))


c(                          map_file,
                          ldetect_file,
                          kg_dir,
                          input_db_f)[!file.exists(c(                          map_file,
                          ldetect_file,
                          kg_dir,
                          input_db_f))]


stopifnot(all(file.exists(c(
                          map_file,
                          ldetect_file,
                          kg_dir,
                          input_db_f))))
