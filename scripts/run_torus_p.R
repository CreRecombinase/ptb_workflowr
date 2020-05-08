library(daprcpp)
library(dplyr)
library(purrr)
library(readr)
library(ldmap)
library(fs)
library(tidyr)
library(stringr)


gf <- snakemake@input[["gwasf"]]
af <- snakemake@input[["annof"]]

torus_path <- snakemake@params[["torus_cmd"]]
stopifnot(!is.null(torus_path))
if (is.null(af)) {
    af <- tempfile()
    write_tsv(tibble::tibble(SNP = character()), af)
}

prior_rf <- snakemake@input[["prior_r"]]

prior_r <- scan(prior_rf, what = character())
od <- snakemake@output[["outputd"]]

run_torus_cmd <- function(gf,af,torus_p=character(0),l1=NA_real_,l2=NA_real_,torus_path=system("which torus")){


  stopifnot(file.exists(gf),
            file.exists(af))
  fo <- fs::file_info(torus_path)
  stopifnot((fo$permissions & "u+x") == "u+x")
  torus_d <- fs::file_temp()
  lik_file <- fs::file_temp()
  if(length(torus_p)>0){
    p_f <- fs::path(torus_d,torus_p,ext="prior")
    stopifnot(!fs::dir_exists(torus_d))
    res_args <- c(
      "-d",
      fs::path_expand(gf),
      "-annot",
      fs::path_expand(af),
      "--load_zval",
      "-lik",
      lik_file,
      "-dump_prior",
      torus_d)
  } else{
    res_args <- c(
      "-d",
      fs::path_expand(gf),
      "-annot",
      fs::path_expand(af),
      "--load_zval",
      "-lik",
      lik_file
    )
  }
  if(!is.na(l1)){
    res_args <- c(res_args,"-l1_lambda",l1)
  }
  if(!is.na(l2)){
    res_args <- c(res_args,"-l2_lambda",l2)
  }
  res <- processx::run(torus_path,args = res_args,echo_cmd = TRUE,echo = TRUE)
  df <- read.table(file = textConnection(res$stdout),skip=1,header=F,stringsAsFactors = F)
  colnames(df) <- c("term", "estimate", "low", "high")

  df <- dplyr::mutate(df,term=stringr::str_replace(term,pattern = "\\.[0-9]+$",replacement = ""),
                      sd=(low-estimate)/(-1.96),z=estimate/sd,p=pnorm(abs(z),lower.tail = FALSE))
  lik <- scan(lik_file,what=numeric())
  file.remove(lik_file)
  df <- tidyr::nest(df, data = tidyr::everything()) %>% dplyr::mutate(lik=lik)
    if( length(torus_p) > 0){
        stopifnot(all(fs::file_exists(p_f)))
        prior_l <- purrr::map(torus_p,function(x){
            fp <- as.character(fs::path(torus_d,x,ext="prior"))
            suppressMessages(
                ret <- vroom::vroom(file = fp,delim = "  ",trim_ws = T,col_names = c("SNP","prior"),col_types = cols("SNP"="c","prior"="d")) %>% dplyr::mutate(region_id=x)
            )
            return(ret)
        })
      fs::file_delete(p_f)
      names(prior_l) <- torus_p
      ret <- list(df=df,priors=prior_l)
  }else{
      ret <- list(df=df)
  }
  return(ret)
}



torus_ret <- run_torus_cmd(gf = gf, af = af, torus_p = prior_r,torus_path=torus_path)

saveRDS(torus_ret$df, snakemake@output[["outputf"]])
if (!dir.exists(od)) {
    fs::dir_create(od, recurse = TRUE)
}
iwalk(torus_ret$priors, function(pr, region_id) {
    trid <- region_id
    pr %>% rename(snp_struct=SNP) %>% mutate(snp_struct=as_ldmap_snp(snp_struct)) %>%
        saveRDS(fs::path(od, trid, ext = "RDS"))
})
