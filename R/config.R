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
    options(clustermq.scheduler = "multicore")
    config_path <-"config/workflow_params_gardner.yaml"
}

if (grepl(x=nodename,pattern =  "dellxps")) {
    config_path <- "~/Dropbox/Repos/ptb_workflowr/config/workflow_params_xps.yaml"
}

allpack <- yaml::read_yaml("config/packages.yaml",handlers=list(package=package_fun))
data_config <- yaml::read_yaml(config_path,
                         handlers=handler_l)
all_feat <- yaml::read_yaml("config/features.yaml")

cload <- purrr::partial(loadd,cache=data_config$cache,envir=parent.frame())
cread <- purrr::partial(readd,cache=data_config$cache)
