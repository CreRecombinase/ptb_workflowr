source("R/handlers.R")

##' read a bed annotation and give it a name
read_annot <- function(feat_file, feat_name){
readr::read_tsv(feat_file, col_names = c("chr", "start", "end"),
             col_types = readr::cols_only(
                 chr = col_character(),
                 start = col_integer(),
                 end = col_integer())) %>%
    dplyr::mutate(name = feat_name)
}



##' Convert a gwas dataframe to a genomicranges
make_range <- function(input_df){
    dplyr::select(input_df, chrom, start = pos) %>%
        mutate(chrom = paste0("chr", chrom),
               start = as.integer(start),
               end = start + 1L) %>%
        GenomicRanges::makeGRangesFromDataFrame()
}



susiefun <- function(bhat,shat,sample_size,L=1,R=diag(length(bhat)),prior=NULL,scaled_prior_variance = 0.1){
susie_bhat(bhat =bhat,
           shat = shat,
           R=R,
           n =sample_size,prior_weights=prior,
           L = L,
             estimate_residual_variance = TRUE,
                                      estimate_prior_variance = FALSE)
}




shim_susie <- function(df,R=diag(nrow(df)),L=1,h_p){
  susiefun(bhat = df$beta,
           shat = df$se,
           sample_size = max(df$N),
           L = L,
           prior=df$prior,
           scaled_prior_variance = h_p*nrow(df)
           )
}

snp_anno_range <- function(snp_anno){
  snp_filename <- fs::path(data_config$anno_dir,snp_anno,ext = "gz")
  trm_field <- rlang::sym(scan(snp_filename,what=character(),sep = "\t",n=2)[2])
  readr::read_tsv(snp_filename) %>%
    filter(!! trm_field != 0L) %>%
    tidyr::separate(SNP,into=c("chrom","pos","ref","alt"),sep=":",remove=F,convert=T) %>%
    mutate(chr=paste0("chr",chrom),start=as.integer(pos),end=start+1L) %>%
    select(chr,start,end) %>%
    GenomicRanges::makeGRangesFromDataFrame()

}


##' Read a bed file
read_anno_r <- function(anno_name,dcf) {
  anno_file <- dcf[str_replace(path_file(dcf),".bed.*","")==anno_name]
  stopifnot(length(anno_file)==1)
  stopifnot(file.exists(anno_file))
  ret <- vroom::vroom(file_in(anno_file), col_names = c("chr", "start", "end"),
               col_types = cols_only(
                 chr = col_character(),
                 start = col_integer(),
                 end = col_integer())) %>%
    GenomicRanges::makeGRangesFromDataFrame()
  attr(ret,"feat_name") <- anno_name
  return(ret)
}





anno_overlap_fun <- function(input_range,gr_df,gw_df,name=input_range@feat_name){

  fr <- GenomicRanges::findOverlaps(gr_df,input_range)
  frovec <- fr@from
  tret <- dplyr::slice(gw_df,frovec) %>% dplyr::select(SNP) %>% mutate(feature=name)

  # ftret <- map2(tret,n_name,function(x,y){
  #   x[[y]] <- rep(1L,nrow(x))
  #   return(x)
  # }) %>% reduce(~full_join(.x,.y,by="SNP")) %>% mutate_at(vars(ends_with("_d")),partial(replace_na,replace=0))
  if(nrow(tret)==0){
    attr(tret,"feature") <- name
  }
  return(tret)
}

bind_results <- function(...){
  input <- list(...)
  map_df(input,~.x$est)
}







torusR <-function(gw_df,anno_l,prior=NA_integer_){
  k <-length(anno_l)
  annomat <-matrix(data=0,nrow=nrow(gw_df),ncol = k)
  p <- nrow(gw_df)
  colnames(annomat) <- map_chr(anno_l,~colnames(.x)[2])
  for(i in 1:seq_along(anno_l)){
    annomat[anno_l[[i]][["SNP"]],i] <- 1L
  }
  use_prior <- all(!is.na(prior))
  res <- daprcpp::torus(locus_id = gw_df$region_id,z_hat = gw_df$`z-stat`,anno_mat = annomat,prior = use_prior,names = colnames(annomat))
  if(use_prior){
    stopifnot(all(prior %in% gw_df$region_id))
    res$est <- mutate(res$est,sd=(low-estimate)/(-1.96),z=estimate/sd,p=pnorm(abs(z),lower.tail = FALSE))
    sub_b <- gw_df$region_id %in% prior
    res$prior <- filter(gw_df,sub_b) %>% mutate(prior=res$prior[sub_b])
  }else{
    res$est <- mutate(res$est,sd=(low-estimate)/(-1.96),z=estimate/sd,p=pnorm(abs(z),lower.tail = FALSE))
  }
  return(res)
}


run_torus_R <- function(gw_df,annomat,prior=NA_integer_){
  use_prior <- all(!is.na(prior))
  res <- daprcpp::torus(locus_id = gw_df$region_id,z_hat = gw_df$`z-stat`,anno_mat = annomat,prior = use_prior,names = colnames(annomat))
  if(use_prior){
    stopifnot(all(prior %in% gw_df$region_id))
    res$est <- mutate(res$est,sd=(low-estimate)/(-1.96),z=estimate/sd,p=pnorm(abs(z),lower.tail = FALSE))
    sub_b <- gw_df$region_id %in% prior
    res$prior <- dplyr::filter(gw_df,sub_b) %>% dplyr::mutate(prior=res$prior[sub_b])
  }else{
    res <- mutate(res,sd=(low-estimate)/(-1.96),z=estimate/sd,p=pnorm(abs(z),lower.tail = FALSE))
  }
  return(res)


}



run_torus_Rdf <- function(gw_df,anno_df,prior=NA_integer_,verbose=F,use_glmnet=TRUE){
  use_prior <- all(!is.na(prior))
  res <- daprcpp::torus_df(locus_id = factor(gw_df$region_id),z_hat = gw_df$`z-stat`,anno_df = anno_df,prior = use_prior,do_verbose = verbose,use_glmnet = use_glmnet)
  if(use_prior){
    stopifnot(all(prior %in% gw_df$region_id))
    res$est <- mutate(res$est,sd=(low-estimate)/(-1.96),z=estimate/sd,p=pnorm(abs(z),lower.tail = FALSE))
    sub_b <- gw_df$region_id %in% prior
    res$prior <- dplyr::filter(gw_df,sub_b) %>% dplyr::mutate(prior=res$prior[sub_b])
  }else{
    res$est <- mutate(res$est,sd=(low-estimate)/(-1.96),z=estimate/sd,p=pnorm(abs(z),lower.tail = FALSE))
  }
  return(res)
}



lik_fun <- function(x){
  x$est
}





forward_op_torus <- function(gw_df,anno_df,f_feat,term_list,i,prior=NA_integer_,verbose=F){
  fn <- term_list[i]
  if(fn %in% f_feat){
    return(set_names(-Inf,fn))
  }
  tdf <- filter(anno_df,feature %in% c(term_list[i],f_feat))
  tret_f <- safely(run_torus_Rdf,otherwise = list(est=list(lik=-Inf)))
  return(set_names(lik_fun(tret_f(gw_df,tdf,verbose = verbose)$result),fn))
}

pr_torus <- function(gw_df,anno_df,feat_v,prior=integer(0)){
  stopifnot(length(prior)>0)
  tdf <- filter(anno_df,feature %in% feat_v)
  retl <-run_torus_Rdf(gw_df,tdf,prior = prior)
  return(list(prior=split(retl$prior,retl$prior$region_id),est=retl$est))
}

combo_lik_df <- function(...){
  bind_rows(list(...)) %>% filter(term!="Intercept")
}

forward_reduce <- function(f_feat,term_list,lik_df){

  stopifnot(nrow(lik_df)==length(term_list))
  bterm <- filter(lik_df,p==min(p)) %>% slice(1)

  stopifnot(!bterm %in% f_feat)
  return(c(f_feat,bterm))
}




forward_reg_torus <- function(res_df,gw_df,anno_df,prior=NA_integer_,iternum=3,predicate=NA_character_,p_thresh=0.05,verbose=F){
  rres_df <- filter(res_df,term!="Intercept",p<p_thresh) %>% arrange(p)
  stopifnot(nrow(rres_df)>0)
  if(!is.na(predicate)){
    stopifnot(is.character(predicate),length(predicate)>0,all(!is.na(predicate)))
    rres_df <-filter(rres_df,str_detect(term,predicate)|str_detect(term,"seq",negate=T))
  }
  term_list <- rres_df$term

  f_feat <-term_list[1]
  term_list <- term_list[-1]

  for( i in seq_len(iternum)){

    if(i==iternum){
      pr <- prior
    }else{
      pr <- NA_integer_
    }

    full_feat <-map(term_list,function(x,y,tb_df,gw_df,pri,verb){
      tdf <- filter(tb_df,feature %in% c(x,y))
      run_torus_Rdf(gw_df,tdf,prior = pri,verbose=verb)
    },y=f_feat,tb_df=anno_df,gw_df=gw_df,pri=pr,verb=verbose)

    all_lik <- map_dbl(full_feat,.f =~lik_fun(.x))
    bf <- which.max(all_lik)
    f_feat <- c(f_feat,term_list[bf])
    best_model <- full_feat[[bf]]
    term_list <- term_list[-bf]
    if(length(term_list)==0){
        return(best_model)
    }
  }
  return(best_model)
}




torusR_df <-function(gw_df,anno_df,prior=NA_integer_){
  return(run_torus_R(gw_df,anno_m,prior=prior))
}



run_torus_p <- function(gwas_filename, anno_filename,torus_d,torus_p) {
  stopifnot(file.exists(gwas_filename))
  stopifnot(file.exists(data_config$torus_path))
  if(length(torus_p)>0){
    stopifnot(!fs::dir_exists(torus_d))
    res_args <- c(
      "-d",
      fs::path_expand(gwas_filename),
      "-annot",
      fs::path_expand(anno_filename),
      "--load_zval",
      "-est",
      "-qtl",
      "-dump_prior",
      torus_d,
      "-regions",
      paste(torus_p,collapse=",")
    )
    p_f <-fs::path(torus_d,torus_p,ext="prior.zstd")
  }else{
    res_args <- c(
      "-d",
      fs::path_expand(gwas_filename),
      "-annot",
      fs::path_expand(anno_filename),
      "--load_zval",
      "-est",
      "-qtl"
    )
  }

  res <- processx::run(data_config$torus_path,args = res_args,echo_cmd = TRUE)
  df <- read.table(file = textConnection(res$stdout),header=T,sep="\t",stringsAsFactors = F)
  colnames(df) <- c("term", "estimate", "low", "high")

}

do_torus_p <- function(feat_name,gr_df,gw_df,gwas_file,annof=as.character(fs::file_temp(ext=".tsv.zstd")),regions=unique(gw_df$region_id)){
  stopifnot(file.exists(gwas_file))
  ur_id=unique(gw_df$region_id)
  stopifnot(all(regions %in% ur_id))
  tf <- archive::file_write(annof,"zstd")
  anno_overlap_fun(input_range = map(feat_name,read_anno_r),
                   gr_df = gr_df,
                   gw_df =gw_df,
                   name = feat_name) %>%
    readr::write_tsv(path = tf)
  td <- fs::path_temp(feat_name)
  retl <- run_torus_p(gwas_filename = gwas_file,
                      anno_filename = annof,
                      torus_d = td,
                      torus_p = regions)
  fs::dir_delete(td)
  file.remove(annof)
  return(retl)
}



do_torus <- function(feat_name,gr_df,gw_df,gwas_file){
  tfn <- fs::file_temp(ext=".tsv.zstd")
  tf <- archive::file_write(tfn,"zstd")
  anno_overlap_fun(input_range = map(feat_name,read_anno_r),
                   gr_df = gr_df,
                   gw_df =gw_df,
                   name = feat_name) %>% readr::write_tsv(tf)

  ret_df <- run_torus(gwas_file,tfn)
  file.remove(tf)
  return(ret_df)
}


gwas_range <- function() {
    db_df <- input_db_f
    snp_df <- dplyr::tbl(dplyr::src_sqlite(path = db_df, create = F), "gwas") %>%
        dplyr::select(SNP = id,
                      chrom = starts_with("ch"),
                      pos) %>%
        mutate(chrom = as.integer(chrom), pos = as.integer(pos)) %>% filter(chrom > 0,chrom < 23)
    ret_df <- group_by(snp_df,chrom) %>% summarise(min_p = min(pos),
                                                   max_p = max(pos)) %>% collect()

}


calc_p <- function(db_df,table_name="gwas"){
  dbc <- dbConnect(drv = MonetDBLite::MonetDBLite(),db_df,create=F)
  db <- src_sql("monetdb",dbc)
  p <- dplyr::tbl(db, table_name)%>% summarise(p=n()) %>% collect() %>% pull(p)
  dbDisconnect(dbc)
  return(p)

}




full_gwas_df<-function(db_df,beta_v="beta", se_v="se", N_v="n",p_v="p",keep_bh_se=TRUE,keep_allele=TRUE,nlines=-1) {
    #dbc <- dbConnect(drv = MonetDBLite::MonetDBLite(),db_df,create=F)
    db <- src_monetdblite(dbdir = db_df,create=F)

    if(nlines>0){
      snp_df <- dplyr::tbl(db, dplyr::sql(glue::glue("SELECT * FROM gwas AS s SAMPLE {nlines}")))
    }else{
      snp_df <-dplyr::tbl(db, "gwas")
    }
    snp_df <- snp_df %>%  dplyr::select(chrom = starts_with("ch"),
                                        pos,
                                        N = !!N_v,
                                        beta = !!beta_v,
                                        se = !!se_v,
                                        p = !!p_v,
                                        a1,
                                        a2
                                        ) %>%
      dplyr::mutate(`z-stat` =  beta/se) %>%
      filter(chrom > 0,chrom < 23)
    if(!keep_bh_se){
      snp_df <- snp_df %>%
        dplyr::select(-beta, -se)
    }
    if(!keep_allele){
      snp_df <- snp_df %>%
        dplyr::select(-a1, -a2)
    }
    snp_df <- snp_df %>%
      dplyr::collect() %>%
      dplyr::distinct(chrom, pos, .keep_all = T) %>%
      dplyr::arrange(chrom, pos)
    dbDisconnect(db$con,shutdown=T)
    return(dplyr::mutate(snp_df,SNP=1:dplyr::n()))
}






assign_reg_df <- function(snp_df,ld_df,max_snp=-1L,min_snp=1L) {
  reg_id <- assign_region(break_chr = ld_df$chrom,
                          break_start = ld_df$start,
                          break_stop = ld_df$stop,
                          break_id = ld_df$region_id,
                          snp_chr = snp_df$chrom,
                          max_size=max_snp,
                          min_size=min_snp,
                          snp_pos = snp_df$pos,
                          assign_all = T)
  dplyr::mutate(snp_df,
                region_id = reg_id)
}




merge_snp_f <- function(geno_f,gwas_df){
  stopifnot(file.exists(geno_f))
  geno_d = snp_attach(geno_f)
  # info_snp <- vroom::vroom(fs::path_ext_set(geno_f,"bim"),delim = "\t",col_names = c("chr",
  #                                                                                    "id",
  #                                                                                    "gpos",
  #                                                                                    "pos",
  #                                                                                    "a0",
  #                                                                                    "a1"
  #                                                                                   ),
  #                          col_types = cols(chr="i",id="c",gpos="i",pos="i","a0"="c",a1="c"))

  info_snp <- dplyr::select(geno_d$map,
                            chr=chromosome,
                            id=marker.ID,
                            pos=physical.pos,
                            a0=allele1,
                            a1=allele2)
  sumstats <- dplyr::select(gwas_df,
                            chr=chrom,
                            pos,
                            ta0=a1,
                            ta1=a2,
                            beta,
                            p
                            ) %>% dplyr::rename(a0=ta0,a1=ta1)

  ret_snp <- dplyr::inner_join(snp_match(sumstats = sumstats,info_snp = info_snp),
                               dplyr::mutate(info_snp,ld_id=1:dplyr::n()))
  dplyr::mutate(dplyr::inner_join(gwas_df,ret_snp,
                        c("chrom"="chr",
                          "pos",
                          "beta",
                          "a1"="a0",
                          "a2"="a1",
                          "p")),SNP=1:dplyr::n())

}

read_ptb_db <- function(db_df, beta_v="beta", se_v="se", N_v="N"){
      dplyr::tbl(dplyr::src_sqlite(path = db_df, create = F), "gwas") %>%
        dplyr::select(SNP = id,
                      chrom = starts_with("ch"),
                      pos,
                      MAJOR = A1,
                      MINOR = A2,
                      N = !!N_v,
                      beta = !!beta_v,
                      se = !!se_v) %>%
    mutate(beta=as.numeric(beta),
           se=as.numeric(se),N=as.numeric(N),
           chrom = as.integer(chrom), pos = as.integer(pos))
}



