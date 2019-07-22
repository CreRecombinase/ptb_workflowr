root_fun <- function(x){
  here::here()
}

path_fun <- function(x){fs::path_expand(do.call(fs::path,x))}


ldd_fun <- function(x){
  readr::read_tsv(x$path,col_types = c(chrom="i",start="i",stop="i",region_id="i"))
}


file_fun <- function(x){
  stopifnot(!is.null(x$path),
            length(x$path)==1)
  if(!is.null(x$check)){
    if(x$check)
      stopifnot(file.exists(x$path))
  }
  return(x$path)
}

files_fun <- function(x){
  stopifnot(!is.null(x$path),
            length(x$path)==1)
  ret <- fs::dir_ls(path=x$path,type = (x$type %||% "any"), glob=x$glob,recurse = (x$recurse %||% F))
  stopifnot(length(ret)>0)
  return(ret)
}
cache_fun <- function(x){
  if(!fs::dir_exists(x$file)){
    fs::dir_create(x$file,recurse = T)
    new_cache(x$file)
  }else{
    drake_cache(x$file)
  }

}

handler_l <- list(
  root = root_fun,
  path = path_fun,
  file = file_fun,
  files = files_fun,
  cache = cache_fun,
  ldd =ldd_fun
)

`%||%` <- function (x, y)
{
    if (is.null(x))
        y
    else x
}

package_fun <- function(x) {

    nx <- x$name
    stopifnot(length(nx)==1)
    if (!(x$load %||% TRUE)) {
        if (!requireNamespace(nx,quietly = TRUE)) {
            if (!is.null(x$github)) {
                if (!require(devtools)) {
                    install.packages("devtools")
                }
                devtools::install_github(x$ghub,ref = x$ref %||%  "master")
            }else if (!is.null(x$bioc)) {
                if (!requireNamespace("BiocManager", quietly = TRUE))
                    install.packages("BiocManager")

                BiocManager::install("GenomicRanges")
            } else {
                install.packages(nx)
            }

        }
      return(requireNamespace(nx,quietly =  TRUE))
    }
    if (suppressWarnings(suppressMessages(!require(nx,character.only = T,quietly = TRUE)))) {
        if (is.null(x$ghub)) {
            install.packages(nx)
            library(nx,character.only = T,quietly =  TRUE)

        }
        if (!require(devtools)) {
            install.packages("devtools")
        }
        devtools::install_github(x$ghub,ref = x$ref %||%  "master")
        return(require(nx,character.only = T,quietly = TRUE))
    }
    return(TRUE)
}
