library(tidyverse)
library(susieR)
library(drake)
library(glue)
library(RSQLite)
library(archive)
library(jsonlite)
library(ldshrink)
library(daprcpp)
library(ldmap)
library(rvest)
library(RSSp)

saveCRDS <- function(object, filename, filter=NULL) {
  stopifnot(filter %in% c('zstd', 'lz4'))
  con = archive::file_write(file = filename, filter=filter)
  open(con)
  saveRDS(object, con)
  close(con)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Read compressed RDS streams
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
readCRDS <- function(filename) {
  con <- archive::file_read(file = filename)
  res <- readRDS(con)
  close(con)
  res
}
