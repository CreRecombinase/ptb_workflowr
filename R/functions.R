source(fs::path(here::here(),"R/handlers.R"))

##' read a bed annotation and give it a name
read_annot <- function(feat_file, feat_name){
readr::read_tsv(feat_file, col_names = c("chr", "start", "end"),
             col_types = readr::cols_only(
                 chr = col_character(),
                 start = col_integer(),
                 end = col_integer())) %>%
    dplyr::mutate(name = feat_name)
}


region_plot <- function(chrom,start,stop,ensembl=biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl"),...){
res <-   biomaRt::getBM(attributes = c('hgnc_symbol',"start_position","end_position"),
      filters = c('chromosome_name','start','end'),
      values = list(chrom,max(start,1L),stop),
      mart = ensembl)
}


match_df <- function(df){
  col_l <- list(character = col_character(),
                list = col_character(),
                double = col_double(),
                integer = col_integer())

  purrr::lift_dl(cols)(purrr::map(df,~col_l[[typeof(.x)]]))
}

prep_sql <- function(url){
  sql_df <- tibble::tibble(lines=readr::read_lines(url))
  sqll <-dplyr::filter(sql_df,
                          str_detect(lines,"^/\\*",negate = TRUE),
                          str_detect(lines,"^--",negate = TRUE),
                          lines!="",
                          str_detect(lines,"DROP TABLE",negate = TRUE),
                          str_detect(lines, "KEY",negate = TRUE)) %>% dplyr::pull(lines) %>% paste0(collapse="")
  sqll <- str_replace(sqll,"ENGINE.+","")
  sqll <- str_replace_all(sqll,"unsigned","")
  str_replace(str_replace(sqll,",\\s*\\)",")"),"` \\(","`(")
}

read_df_ucsc <- function(tblname,
                         sql_url=glue::glue("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/{tblname}.sql"),
                         df_url=glue::glue("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/{tblname}.txt.gz")){
  sql_st <- prep_sql(sql_url)
  dbf <- fs::file_temp()
  db <- src_sqlite(dbf,create=TRUE)
  DBI::dbSendQuery(db$con,sql_st)
  kgdb <- tbl(db,db_list_tables(db$con)[1])
  db_df <- readr::read_tsv(df_url,col_names=colnames(kgdb),col_types = match_df(collect(kgdb)))
  DBI::dbDisconnect(db$con)
  fs::file_delete(dbf)
  return(db_df)
}


ucsc2monetdb <- function(tblname,dbdir=fs::file_temp()){
  db <- MonetDBLite::src_monetdblite(dbdir=dbdir)
  dbWriteTable(db$con,tblname,read_df_ucsc(tblname))
  dbDisconnect(db$con)
  return(dbdir)
}



gene2genebody <- function(df){

  df$start_str <-  map(str_split(str_replace(df$exonStarts,",$",""),","),as.integer)
  df$stop_str <- map(str_split(str_replace(df$exonEnds,",$",""),","),as.integer)

  # stopifnot(all.equal(lengths(start_str),lengths(stop_str)),all.equal(lengths(start_str),df$exonCount))

  df <- dplyr::select(df,-exonStarts,-exonEnds,-exonCount)
  df$exon_l <- pmap(df,function(start_str,stop_str,txStart,txEnd,cdsStart,cdsEnd,...){
      if(cdsStart>txStart){
        pref_tbl <- tibble::tibble(start=txStart,end=cdsStart,type="UTR")
      }
      else{
        pref_tbl <- tibble::tibble(start=integer(),end=integer(),type=character())
      }
      if(cdsEnd<txEnd){
        post_tbl <- tibble::tibble(start=cdsEnd,end=txEnd,type="UTR")
      }
      else{
        post_tbl <- tibble::tibble(start=integer(),end=integer(),type=character())
      }
      dplyr::bind_rows(pref_tbl,tibble::tibble(start=start_str,end=stop_str,type="CDS"),post_tbl )

    })
  unnest(df,exon_l) %>% arrange(as.integer(str_replace(chrom,"chr","")),start,end)
}


nearest_gene_df <-function(snpdf,dbdir,needs_mRNA=TRUE){
  osnpdf <- mutate(snpdf,width=1L,seqnames=glue("chr{chrom}"),strand = "*",start=pos) %>%plyranges::as_granges()
  MonetDBLite::monetdblite_shutdown()
  db <- MonetDBLite::src_monetdblite(dbdir=dbdir)
  ref <- tbl(db,"kgXref")
  if(needs_mRNA){
    ref <- dplyr::filter(ref,!is.na(mRNA))
  }
  kgdf <- ref %>% dplyr::select(kgID,geneSymbol,description) %>% distinct() %>%
    inner_join(tbl(db,"knownGene"),by=c("kgID"="name")) %>%
    dplyr::mutate(width=txEnd-txStart) %>%
    dplyr::rename(start=txStart,seqnames=chrom) %>%
    plyranges::as_granges()
  ret_df <- left_join(snpdf,plyranges::join_nearest(osnpdf, kgdf) %>%
    as_tibble() %>%
  dplyr::select(rsid,geneSymbol,description))
  dbDisconnect(db$con)
  return(ret_df)
}


gene_body_df <- function(chr,start,end,dbdir,needs_mRNA=TRUE){
  MonetDBLite::monetdblite_shutdown()
  db <- MonetDBLite::src_monetdblite(dbdir=dbdir)
  if(typeof(chr)=="integer")
    chr <- paste0("chr",chr)
  stopifnot(substr(chr,1,3)=="chr",length(chr)==1)
  ref <- tbl(db,"kgXref")
  if(needs_mRNA){
    ref <- dplyr::filter(ref,!is.na(mRNA))
  }
  kgdf <- tbl(db,"knownGene") %>%
    dplyr::filter(chrom==chr,!(txEnd<start), !(txStart>end)) %>%
    dplyr::inner_join(ref,by=c("name"="kgID")) %>%
    dplyr::collect()
  if(nrow(kgdf)==0){
    dbDisconnect(db$con)
    return(tibble::tibble(chrom=character(),
                          strand=character(),
                          txStart=integer(),
                          txEnd=integer(),
                          cdsStart=integer(),
                          cdsEnd=integer(),
                          geneSymbol=character(),
                          description=character(),
                          start=integer(),
                          end=integer(),
                          type=character()))
  }
  kg <-   kgdf %>%  gene2genebody()  %>%
    dplyr::select(chrom,
                  strand,
                  txStart,
                  txEnd,
                  cdsStart,
                  cdsEnd,
                  geneSymbol,
                  description,
                  start,
                  end,
                  type)
  dbDisconnect(db$con)
  return(kg)
}


##' Convert a gwas dataframe to a genomicranges
make_range <- function(input_df){
    dplyr::select(input_df, chrom, start = pos,dplyr::everything()) %>%
        mutate(chrom = paste0("chr", chrom),
               start = as.integer(start),
               end = start + 1L) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
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




shim_rssp <- function(df,L=4,h_p,geno_f=NULL){
  if(!is.null(geno_f)){
    stopifnot(file.exists(geno_f))
    geno_d <- snp_attach(geno_f)
    R <- cor(geno_d$genotypes[,df$ld_id])
    evdR <- eigen(R)
    Q <- evdR$vectors
    D <- evdR$values
  }else{
    Q <-diag(nrow(df))
    D <- rep(1.0,nrow(df))
  }
  quh <- RSSp::convert_quh(uhat = df$beta/df$se,Q = Q)
  RSSp::RSSp_estimate(quh=quh,D=D,sample_size=max(df$N),trait_id = df$region_id[1])
}

shim_susie <- function(df,L=4, h_p, geno_f=NULL){
  if(is.null(df$p)){
    df$p <- 2*pnorm(abs(df$beta/df$se),lower.tail=F)
  }
  top_id <- which.min(df$p)
  if(!is.null(geno_f)){
    stopifnot(file.exists(geno_f))
    geno_d <- snp_attach(geno_f)
    R <- cor(geno_d$genotypes[,df$ld_id])
  }else{
    R <- Matrix::Diagonal(nrow(df))
  }
  Rd <- R[,top_id]
  res <- susiefun(bhat = df$beta,
           shat = df$se,
           R = R,
           sample_size = max(df$N),
           L = L,
           prior = df$prior,
           scaled_prior_variance = h_p*nrow(df)
  )
  return(list(df=dplyr::mutate(df,t_r=Rd,pip=res$pip),
              susie_res=res))
}


merge_susie <- function(df,susie_res){
  dplyr::mutate(df,pip=susie_res$pip)
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
read_anno_r <- function(anno_file,anno_name) {
#  anno_file <- dcf[str_replace(path_file(dcf),".bed.*","")==anno_name]
  stopifnot(length(anno_file) == 1)
  stopifnot(file.exists(anno_file))
  ret <- vroom::vroom(anno_file, col_names = c("chr", "start", "end"),
               col_types = cols_only(
                 chr = col_character(),
                 start = col_integer(),
                 end = col_integer())) %>%
    GenomicRanges::makeGRangesFromDataFrame()
  attr(ret,"feat_name") <- anno_name
  return(ret)
}





anno_overlap_fun <- function(input_range,gr_df,name=input_range@feat_name){

    fr <- GenomicRanges::findOverlaps(gr_df,input_range)
    sub_df <- IRanges::subsetByOverlaps(gr_df,input_range)
    tret <- tibble::tibble(SNP = GenomicRanges::mcols(sub_df)[["SNP"]],feature = name)
#    frovec <- fr@from

 # tret <- dplyr::slice(gw_df,frovec) %>% dplyr::select(SNP) %>% mutate(feature=name)

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




full_gwas_df<-function(db_df,beta_v="beta", se_v="se", N_v="N",p_v="pval",keep_bh_se=TRUE,keep_allele=TRUE,nlines=-1) {
                                        #dbc <- dbConnect(drv = MonetDBLite::MonetDBLite(),db_df,create=F)
    my_c <- cols(
        id = col_skip(),
        chr = col_integer(),
        pos = col_integer(),
        A1 = col_character(),
        A2 = col_character(),
        N = col_double(),
        freq = col_skip(),
        beta = col_skip(),
        se = col_skip(),
        pval = col_skip(),
        Q = col_skip(),
        het = col_skip(),
        N.local = col_skip(),
        freq.local = col_skip(),
        beta.local = col_skip(),
        se.local = col_skip(),
        pval.local = col_skip(),
        N.23andMe = col_skip(),
        freq.23andMe = col_skip(),
        beta.23andMe = col_skip(),
        se.23andMe = col_skip(),
        pval.23andMe = col_skip()
    )

    for (cn in c(beta_v,se_v,N_v,p_v)) {
        my_c[["cols"]][[cn]] <- col_double()
    }
    snp_df <- vroom::vroom(db_df,delim="\t",col_types = my_c) %>%  dplyr::select(chrom = chr,
                                                                                 pos,
                                                                                 N = !!N_v,
                                                                                 beta = !!beta_v,
                                                                                 se = !!se_v,
                                                                                 p = !!p_v,
                                                                                 a1 = A1,
                                                                                 a2 = A2)
    if(nlines>0){
      snp_df <- dplyr::sample_n(nlines,replace = F)
    }else{

    }
    snp_df <- snp_df  %>%
      dplyr::mutate(`z-stat` =  beta/se) %>%
      dplyr::filter(chrom > 0,chrom < 23)
    if(!keep_bh_se){
      snp_df <- snp_df %>%
        dplyr::select(-beta, -se)
    }
    if(!keep_allele){
      snp_df <- snp_df %>%
        dplyr::select(-a1, -a2)
    }
    snp_df <- snp_df %>%
      dplyr::distinct(chrom, pos, .keep_all = T) %>%
      dplyr::arrange(chrom, pos)

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




merge_snp_f <- function(geno_f,gwas_df,strand_flip=TRUE){
  stopifnot(file.exists(geno_f))
  geno_d = snp_attach(geno_f)


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
                            ) %>%
    dplyr::rename(a0=ta0,a1=ta1)

  ret_snp <- dplyr::inner_join(snp_match(sumstats = sumstats,info_snp = info_snp,strand_flip=strand_flip),
                               dplyr::mutate(info_snp,ld_id=1:dplyr::n()))
  dplyr::select(gwas_df,-beta,-`z-stat`,-a1,-a2,-p) %>%
    dplyr::inner_join(ret_snp,
                      c("chrom"="chr",
                        "pos")) %>%
    mutate(SNP=1:dplyr::n(),`z-stat`=beta / se)
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

write_gwas <- function(gwas_df,gf=tempfile(fileext=".txt.gz")){
  fs::dir_create(fs::path_dir(gf))
  dplyr::select(gwas_df,SNP,region_id,`z-stat`) %>% write_tsv(path=gf)
  return(gf)
}


write_anno <- function(anno_df=tibble(SNP=integer(),feature=character()),p=max(max(anno_df$SNP),1L),af=tempfile(fileext=".txt.gz")){
  if(is.null(anno_df)){
    return(write_anno(af=af))
  }
  fs::dir_create(fs::path_dir(af))
  spread_anno_l <- make_matrix(p = p,anno_df = anno_df)
  spread_anno_df <- tibble::as_tibble(magrittr::set_colnames(spread_anno_l$annomat,paste0(spread_anno_l$names,"_d"))) %>%
    dplyr::mutate(SNP=1:dplyr::n()) %>%
    dplyr::select(SNP,dplyr::everything()) %>%
    dplyr::filter_at(.vars = dplyr::vars(-SNP),dplyr::any_vars(. != 0))
  readr::write_tsv(spread_anno_df,path=af)
  return(af)

}






