config_basef <- "../workflow/config_base.yaml"
annotf <- "../workflow/annots.yaml"

host_fun <- function(x) {
  xo <- x$options
  for (n in names(xo)) {
    if (fs::dir_exists(xo[[n]])) {
      return(n)
    }
  }
}

dep_fun <- function(x) {
  if (is.null(x$pref)) 
    return(x$path[[x$host]])
  return(fs::path(x$pref,x$path[[x$host]]))
}

handler_l <- list(
  Host = host_fun,
  Dep = dep_fun
)

config <- yaml::read_yaml(
                    config_basef,
                  handlers =  handler_l)
config_d <- config$paths
annot <- yaml::read_yaml(annotf)
