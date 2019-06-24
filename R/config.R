read_yaml("config/packages.yaml",handlers=list(package=package_fun))

nodename <- Sys.info()["nodename"]
if(str_detect(nodename,"helab")){
  config_path <-"config/workflow_params_desktop.json"
  options(clustermq.scheduler = "multicore")
}
if(str_detect(nodename,"midway2")){
    options(
    clustermq.scheduler = "slurm",
    clustermq.template = "/scratch/midway2/nwknoblauch/ptb_workflowr/slurm_clustermq.tmpl"
)
    config_path <-"config/workflow_params_rcc.yaml"
}

if (str_detect(nodename, "dellxps")) {
    config_path <- "~/Dropbox/Repos/ptb_workflowr/config/workflow_params_xps.yaml"
}


data_config <- read_yaml(config_path,
                         handlers=handler_l)


cload <- purrr::partial(loadd,cache=data_config$cache,envir=parent.frame())
cread <- purrr::partial(readd,cache=data_config$cache)

# h <- curl::new_handle()
# curl::handle_setopt(handle = h,httpauth=2,
#                     userpwd="mod:12apples",verbose=T)
#
# url <- "https://mnlab.uchicago.edu/mod/tmp/peaks/"
# dl_files <- curl::curl_fetch_memory(url,handle = h) %>% magrittr::extract2("content") %>%
# xml2::read_html() %>%
#   html_node("table") %>%html_nodes("a") %>% html_attr("href")
# dl_files <- dl_files[str_detect(dl_files,pattern = ".*bed.bz2")]
#
# pool <- curl::new_pool()
# cb <- function(req){cat("done:", req$url, ": HTTP:", req$status, "\n")}
# destdir <- "/home/nwknoblauch/Dropbox/scratch/ptb_scratch/new_bed/"
# map(dl_files,~curl::curl_download(paste0(url,.x),destfile = fs::path(destdir,.x),handle = h))
#

