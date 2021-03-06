source(fs::path(here::here(),"R","functions.R"))


nodename <- Sys.info()["nodename"]
if(grepl(x = nodename,pattern = "helab")){
  config_path <-c("config/data/desktop.yaml","config/params/desktop.yaml")
  options(clustermq.scheduler = "multicore")
}
if(grepl(x=nodename,pattern="midway2")){
  #     options(
  #     clustermq.scheduler = "slurm",
  #     clustermq.template = "/scratch/midway2/nwknoblauch/ptb_workflowr/slurm_clustermq.tmpl"
  # )
  config_path <-c("config/data/rcc.yaml","config/params/rcc.yaml")
}


if(grepl(x=nodename,pattern="cri")){
  options(unzip = "/usr/bin/unzip")
  options(clustermq.scheduler = "multicore")
  config_path <-"config/workflow_params_gardner.yaml"
}

if (grepl(x = nodename,
          pattern =  "dellxps")) {
    config_path <- c(
        "config/data/xps.yaml",
        "config/params/xps.yaml"
    )
}

allpack <- yaml::read_yaml(fs::path(here::here(),"config/packages.yaml"),
                           handlers = list(package = package_fun))
data_config <- map(config_path, ~yaml::read_yaml(fs::path(here::here(), .x),
                               handlers=handler_l)) %>% flatten()
feat_l <- yaml::read_yaml(fs::path(here::here(),"config/features.yaml"),
                          handlers = list(model = model_fun))
all_feat <- feat_l$features
model_df <- bind_rows(feat_l$models)
all_feat <- unique(c(unlist(model_df$features), all_feat))

gf <- as.character(fs::path(data_config$data$scratch, paste0("gwas_f"),ext = "txt.gz"))
stopifnot(fs::dir_exists(fs::path_dir(gf)),all(nchar(gf) > 0))

af <- fs::path(data_config$data$scratch,"new_bed",all_feat,ext = "bed.bz2")
#af <- unique(paste0(str_replace(af, ".bed$", ""), ".txt.gz"))
#af <- af[str_replace(fs::path_file(af), ".txt.gz", "") %in% all_feat]

## out_f  <- cross2(models,all_feat) %>%  map_chr(~fs::path(fs::file_temp(), .x[[1]], .x[[2]],ext = ".txt.gz"))


## cload <- purrr::partial(loadd,
##                         cache = data_config$cache,
##                         envir = parent.frame())
## cread <- purrr::partial(readd,
##                         cache = data_config$cache)
3147274556
