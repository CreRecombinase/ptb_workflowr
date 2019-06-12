##' read a bed annotation and give it a name
read_annot <- function(feat_file, feat_name){
readr::read_tsv(feat_file, col_names = c("chr", "start", "end"),
             col_types = readr::cols_only(
                 chr = col_character(),
                 start = col_integer(),
                 end = col_integer())) %>%
    dplyr::mutate(name = feat_name)
}

                                        #https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE85330&format=file

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


read_scz <- function(){
snp_df <- readr::read_tsv(file_in("~/Downloads/SCZ/Summary_statistics.gz"),col_types = cols_only(chr = col_character(),
  bp = col_integer(),
  z = col_double())) %>%
  rename(pos=bp) %>%  mutate(chrom=as.integer(gsub("chr","",chr)))
  reg_id <- assign_region(break_chr = ld_df$chrom,
                            break_start = ld_df$start,
                            break_stop = ld_df$stop,
                            break_id = ld_df$region_id,
                            snp_chr = snp_df$chrom,
                            snp_pos = snp_df$pos,
                            assign_all = T)
    dplyr::mutate(snp_df,
                  region_id = reg_id) %>%
      tidyr::unite(SNP,chrom,pos,sep="_",remove=F) %>%
      dplyr::select(SNP,chrom,pos,region_id,`z-stat`=z)

}



##' Read a bed file
read_anno_r <- function(anno_name) {
    anno_file <- fs::dir_ls(path = data_config$anno_dir,recursive = T,glob = paste0("*",anno_name,".bed*"))
  stopifnot(length(anno_file)==1)
    vroom::vroom(anno_file, col_names = c("chr", "start", "end"),
                 col_types = cols_only(
                   chr = col_character(),
                   start = col_integer(),
                   end = col_integer())) %>%
      GenomicRanges::makeGRangesFromDataFrame()
}


shuffle_anno <- function(input_range,range_df){

    ilist_df <- GenomicRanges::split(input_range,GenomicRanges::seqnames(input_range)) %>% S4Vectors::sapply(GenomicRanges::width)  %>% enframe(name = "chrom",value = "widths")
    n_range_df <- mutate(range_df,chrom = paste0("chr",chrom)) %>%
        inner_join(ilist_df) %>%
        unnest()
        new_df <- pmap_dfr(n_range_df,function(chrom,min_p,max_p,widths) {
            nmax_p <- max_p - widths
            tibble::tibble(chrom = chrom,start = sample(min_p:nmax_p,1),end = start + widths)
            }) %>% GenomicRanges::makeGRangesFromDataFrame()

}




anno_overlap_fun <- function(input_range,gr_df,gw_df,name){
    n_name <- paste0(name, "_d")
    fr <- purrr::map(input_range,~GenomicRanges::findOverlaps(gr_df,.x))
    frovec <- map(fr,~.x@from)
    tret <- map(frovec,~dplyr::slice(gw_df,.x) %>% dplyr::select(SNP))

    ftret <- map2(tret,n_name,function(x,y){
      x[[y]] <- rep(1L,nrow(x))
      return(x)
    }) %>% reduce(~full_join(.x,.y,by="SNP")) %>% mutate_at(vars(ends_with("_d")),partial(replace_na,replace=0))

    return(ftret)
}

bind_annotations <- function(...){
  input <- list(...)
  map_df(input,~tidyr::gather(.x,key="feature",value="value",-SNP) %>% select(-value))
}



#forward_torus <- function(result_df,anno_df)




torusR <-function(gw_df,anno_l,prior=NA_integer_){
  k <-length(anno_l)
  p <- nrow(gw_df)
  annomat <-matrix(data=0,nrow=p,ncol = k)
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
    res <- mutate(res,sd=(low-estimate)/(-1.96),z=estimate/sd,p=pnorm(abs(z),lower.tail = FALSE))
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
    res$prior <- filter(gw_df,sub_b) %>% mutate(prior=res$prior[sub_b])
  }else{
    res <- mutate(res,sd=(low-estimate)/(-1.96),z=estimate/sd,p=pnorm(abs(z),lower.tail = FALSE))
  }
  return(res)


}

forward_reg_torus <- function(res_df,gw_df,anno_df,prior=NA_integer_,iternum=3){
  rres_df <- filter(res_df,term!="Intercept",p<0.05) %>% arrange(p)
  f_feat <-slice(rres_df,1) %>% pull(term)
  n_feat <-slice(rres_df,2) %>% pull(term)
  new_feat <- filter(anno_df,feature %in% c(f_feat,n_feat))
  rres_df <- slice(rres_df,-c(1,2))
  new_res <-torusR_df(gw_df,new_feat)
  for(i in 1:iternum){
    newr_res <- filter(new_res,term!="Intercept",p<0.05) %>% arrange(p) %>% mutate(term=as.character(term))
    f_feat <- c(newr_res$term)
    n_feat <-slice(rres_df,1) %>% pull(term)
    stopifnot(!n_feat %in% f_feat)
    rres_df <- slice(rres_df,-c(1))
    new_feat_df <- filter(anno_df,feature %in% c(f_feat,n_feat))
    if(i==iternum){
      new_res <-torusR_df(gw_df,new_feat_df,prior = prior)
      new_res$prior_l <- split(new_res$prior,new_res$prior$region_id)
    }else{
      new_res <-torusR_df(gw_df,new_feat_df)
    }
  }
  return(new_res)
}

torusR_df <-function(gw_df,anno_df,prior=NA_integer_){
  anno_df <- mutate(anno_df,feature=factor(feature))
  num_feat <- length(levels(anno_df$feature))
  p <- nrow(gw_df)
  anno_m <- matrix(0,nrow=p,ncol=num_feat)
  colnames(anno_m) <- levels(anno_df$feature)
  anno_m[cbind(anno_df$SNP,as.integer(anno_df$feature))] <- 1L
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
  if(length(torus_p)>0){
    stopifnot(all(fs::file_exists(p_f)))
    prior_l <- map(torus_p,function(x){
      fp <- as.character(fs::path(torus_d,x,ext="prior.zstd"))
      fa <- archive::file_read(fp)
      suppressMessages(
        ret <- read_tsv(file = fa,col_names = c("SNP","prior"),col_types = cols("SNP"="i","prior"="d")) %>% mutate(region_id=x)
      )
      return(ret)
    })
    names(prior_l) <- torus_p
    return(list(df=df,priors=prior_l))
  }else{
    return(list(df=df))
  }
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



full_gwas_df<-function(beta_v, se_v, N_v,keep_bh_se=TRUE) {
    db_df <- input_db_f
    snp_df <- dplyr::tbl(dplyr::src_sqlite(path = db_df, create = F), "gwas") %>%
        dplyr::select(SNP = id,
                      chrom = starts_with("ch"),
                      pos,
                      MAJOR = A1,
                      MINOR = A2,
                      N = !!N_v,
                      beta = !!beta_v,
                      se = !!se_v) %>%
        dplyr::mutate(beta=as.numeric(beta),
                      se=as.numeric(se),
                   `z-stat` =  beta/se,
                   N = as.numeric(N)) %>%
        mutate(chrom = as.integer(chrom), pos = as.integer(pos)) %>% filter(chrom > 0,chrom < 23)
    if(!keep_bh_se){
    snp_df <- snp_df %>%   dplyr::select(-beta, -se)
    }
    snp_df <- snp_df %>%
    dplyr::collect() %>%
      tidyr::unite(col = "allele", MAJOR, MINOR, sep = ",") %>%
      dplyr::distinct(chrom, pos, .keep_all = T) %>%
      arrange(chrom, pos)

    reg_id <- assign_region(break_chr = ld_df$chrom,
                            break_start = ld_df$start,
                            break_stop = ld_df$stop,
                            break_id = ld_df$region_id,
                            snp_chr = snp_df$chrom,
                            snp_pos = snp_df$pos,
                            assign_all = T)
    dplyr::mutate(snp_df,
                  region_id = reg_id)
}






assign_reg_df <- function(snp_df,ld_df) {
    reg_id <- assign_region(break_chr = ld_df$chrom,
                            break_start = ld_df$start,
                            break_stop = ld_df$stop,
                            break_id = ld_df$region_id,
                            snp_chr = snp_df$chrom,
                            snp_pos = snp_df$pos,
                            assign_all = T)
    mutate(snp_df,
           region_id = reg_id)
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



