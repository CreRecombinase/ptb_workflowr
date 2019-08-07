source("R/functions.R")

nodename <- Sys.info()["nodename"]
if(grepl(x = nodename,pattern = "helab")){
  config_path <-"config/workflow_params_desktop.yaml"
  options(clustermq.scheduler = "multicore")
}
if(grepl(x=nodename,pattern="midway2")){
#     options(
#     clustermq.scheduler = "slurm",
#     clustermq.template = "/scratch/midway2/nwknoblauch/ptb_workflowr/slurm_clustermq.tmpl"
# )
    config_path <-"config/workflow_params_rcc.yaml"
}


if(grepl(x=nodename,pattern="cri")){
    options(unzip = "/usr/bin/unzip")
    options(clustermq.scheduler = "multicore")
    config_path <-"config/workflow_params_gardner.yaml"
}

if (grepl(x = nodename,
          pattern =  "dellxps")) {
    config_path <- "~/Dropbox/Repos/ptb_workflowr/config/workflow_params_xps.yaml"
}

allpack <- yaml::read_yaml("config/packages.yaml",
                           handlers = list(package = package_fun))
data_config <- yaml::read_yaml(config_path,
                               handlers=handler_l)
feat_l <- yaml::read_yaml("config/features.yaml",
                          handlers = list(model = model_fun))
all_feat <- feat_l$features
model_df <- bind_rows(feat_l$models)
all_feat <- unique(c(unlist(model_df$features), all_feat))

gf <- as.character(fs::path(data_config$data$home,
                            "ptb_scratch",
                            "gwas_f",
                            ext = "txt.gz"))

af <- fs::path_ext_remove(data_config$data$anno)
af <- paste0(str_replace(af, ".bed$", ""), ".txt.gz")
af <- af[str_replace(fs::path_file(af), ".txt.gz", "")
         %in% all_feat]



cload <- purrr::partial(loadd,
                        cache=data_config$cache,
                        envir=parent.frame())
cread <- purrr::partial(readd,
                        cache=data_config$cache)
