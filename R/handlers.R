root_fun <- function(x){
  here::here()
}

path_fun <- function(x){fs::path_expand(do.call(fs::path,x))}

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
  cache = cache_fun
)

package_fun <- function(x) {
    if (!require(x,character.only = T)) {
        if (is.null(x$ghub)) {
            install.packages(x)
            library(x,character.only = T)
            return(TRUE)
        }
        install_github(x$ghub,ref = x$ref %||%  "master")
        library(x,character.only = T)
    }
}
