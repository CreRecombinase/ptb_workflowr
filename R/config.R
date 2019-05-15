nodename <- Sys.info()["nodename"]
if(str_detect(nodename,"helab")){
  config_path <-"config/workflow_params_desktop.json"
}
if(str_detect(nodename,"midway2")){
    options(
    clustermq.scheduler = "slurm",
    clustermq.template = "/scratch/midway2/nwknoblauch/ptb_workflowr/slurm_clustermq.tmpl"
)
    config_path <-"config/workflow_params_rcc.json"

}

if (str_detect(nodename, "XPS")) {
    config_path <- "config/workflow_params_xps.json"
}

data_config <- jsonlite::read_json(config_path)
