anno_files <- snakemake@input[["annof"]]
anno_names <- snakemake@params[["anno_names"]]
regid <- snakemake@params[["region_id"]]

  library(ldmap)
  library(stringr)
  library(EigenH5)
  library(tidyr)
  library(vroom)
  library(purrr)
  library(dplyr)

  anno_dff <- function(anno_files, anno_names, ldmrid) {
    map2_dfr(anno_files, anno_names, function(x, y) {
      read_bed(x, read_fun = purrr::partial(vroom::vroom, delim = "	")) %>%
        dplyr::filter(region_in_region(ldmap_region, ldetect_EUR,TRUE) == ldmrid) %>%
        mutate(anno = y)
    }) %>% mutate(ldmap_region=as_ldmap_region(ldmap_region))
  }




  stopifnot(!is.null(regid))


  anno_result_df <- anno_dff(anno_files, anno_names, as.integer(regid))



  saveRDS(anno_result_df, snakemake@output[["anno"]])
