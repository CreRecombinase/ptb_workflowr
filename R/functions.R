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
           L = L)
}


shim_susie <- function(df,R=diag(nrow(df)),L=1){
  susiefun(bhat = df$beta,shat = df$se,sample_size = max(df$N),L = L)
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
    anno_file <- fs::path(data_config$anno_dir,anno_name,ext = "bed")
    readr::read_tsv(anno_file, col_names = c("chr", "start", "end"),
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





eqtl_feather <- function(tissue_name){
  eqtl_file <- fs::path(data_config$eqtl_feather_dir,tissue_name,ext = "feather")
  feather::read_feather(eqtl_file,columns = c("chrom","pos")) %>%rename(start=pos) %>%
    mutate(chrom=glue::glue("chr{chrom}"),end=start+1L) %>% GenomicRanges::makeGRangesFromDataFrame()
}

read_meta_qtl <- function(eqtl_path){
  headers <- scan(eqtl_path,what=character(),nlines = 1)
  cols <- headers[toupper(headers)==headers]
  bm_df <- vroom::vroom(eqtl_path,n_max = 10)
}



write_bedf <- function(gr,path){
  df <- data.frame(seqnames=GenomicRanges::seqnames(gr),
  starts=GenomicRanges::start(gr)-1,
  ends=GenomicRanges::end(gr),
  names=c(rep(".", length(gr))),
  scores=c(rep(".",length(gr))),
  strands=GenomicRanges::strand(gr))

  write.table(df, file=path, quote=F, sep="\t", row.names=F, col.names=F)

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

write_snp <- function(gwas_df,filename) {
    dplyr::select(gwas_df,SNP,locus = region_id,`z-stat`) %>%
        readr::write_tsv(gwas_df,drake::file_out(filename))
}

write_anno <- function(anno_df,filename) {
    readr::write_tsv(anno_df,drake::file_out(filename))
}





run_torus <- function(gwas_filename, anno_filename) {
    res_args <- c(
        "-d",
        fs::path_expand(gwas_filename),
        "-annot",
        fs::path_expand(anno_filename),
        "--load_zval",
        "-est",
        "-qtl")
    res <- processx::run(data_config$torus_path,args = res_args,echo_cmd = TRUE)
    res_x <- read.table(textConnection(res$stdout),stringsAsFactors = F)
    fc <- which(res_x$V1=="1")
    stopifnot(length(fc)==1)
    df <- res_x[1:(fc-1), ]
    colnames(df) <- c("term", "estimate", "low", "high")
    return(df)
}


do_torus <- function(feat_name,gr_df,gw_df,gwas_file){
  tf <- fs::file_temp(ext=".tsv.gz")
  anno_overlap_fun(input_range = map(feat_name,read_anno_r),
                   gr_df = gr_df,
                   gw_df =gw_df,
                   name = feat_name) %>% readr::write_tsv(tf)
  ret_df <- run_torus(gwas_file,tf)
  file.remove(tf)
  return(ret_df)
}


do_torus_snp <- function(feat_name,gr_df,gw_df,gwas_file){
  tf <- fs::file_temp(ext=".tsv.gz")
  anno_overlap_fun(input_range = snp_anno_range(feat_name),
                   gr_df = gr_df,
                   gw_df =gw_df,
                   name = feat_name) %>%
    readr::write_tsv(tf)
  ret_df <- run_torus(gwas_file,tf)
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



full_gwas_df<-function(beta_v, se_v, N_v) {
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
        dplyr::mutate(
                   `z-stat` = as.numeric(beta) / as.numeric(se),
                   N = as.numeric(N)) %>%
        mutate(chrom = as.integer(chrom), pos = as.integer(pos)) %>% filter(chrom > 0,chrom < 23)%>%
        dplyr::select(-beta, -se) %>%
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


sub_reg_df  <- function(df,t_chrom, t_start, t_stop){
  df %>%
        dplyr::filter(chrom == t_chrom, between(pos, t_start, t_stop)) %>%
        dplyr::collect() %>%
    tidyr::unite(col = "allele", MAJOR, MINOR, sep = ",") %>%
    dplyr::distinct(chrom, pos, .keep_all = T) %>%
    arrange(chrom, pos)
}



snp_reg <- function(t_chrom, t_start, t_stop, beta_v, se_v, N_v){
    db_df <- input_db_f
    dplyr::tbl(dplyr::src_sqlite(path = db_df, create = F), "gwas") %>%
        dplyr::select(SNP = id,
                      chrom = starts_with("ch"),
                      pos,
                      MAJOR = A1,
                      MINOR = A2,
                      N = !!N_v,
                      beta = !!beta_v,
                      se = !!se_v) %>%
        dplyr::mutate(
                   `z-stat` = as.numeric(beta) / as.numeric(se),
                   N = as.numeric(N)) %>%
        mutate(chrom = as.integer(chrom), pos = as.integer(pos)) %>%
        dplyr::select(-beta, -se) %>%
        dplyr::filter(chrom == t_chrom, between(pos, t_start, t_stop)) %>%
        dplyr::collect() %>%
    tidyr::unite(col = "allele", MAJOR, MINOR, sep = ",") %>%
    dplyr::distinct(chrom, pos, .keep_all = T) %>%
    arrange(chrom, pos)
}

map_snp_reg <- function(df, beta_v, se_v, N_v){
    rename(df, t_chrom = chrom, t_start = start, t_stop = stop) %>%
        select(-region_id) %>%
        pmap(snp_reg, beta_v = beta_v, se_v = se_v, N_v = N_v)
}


read_map <- function(chrom, start, stop, ...){

  snp_df_c <- EigenH5::read_vector_h5(map_file, "SNPinfo/chr")

  which_c <-which(snp_df_c == chrom)
  stopifnot(length(which_c) > 0)
  EigenH5::read_df_h5(map_file, "SNPinfo", subcols = c("pos", "map"), subset = which_c) %>%
      filter(dplyr::between(pos, start, stop)) %>% distinct(map, .keep_all=T) %>%
      arrange(pos)
}

map_read_map <- function(df){
    purrr::pmap(df, read_map)
}





read_df_reg <- function(pop, chrom, start, stop, read_map = F){

    input_f <- fs::path(kg_dir,
                       glue::glue("{pop}.chr{chrom}.h5"))
  snp_df_p <- EigenH5::read_vector_h5(input_f, "SNPinfo/pos")
  subcols <- c("SNP", "allele", "chr", "pos", "snp_id")
  if (read_map){
  subcols <- c(subcols, "map")
  }
  EigenH5::read_df_h5(input_f, "SNPinfo", subcols = subcols,
             subset = which(dplyr::between(snp_df_p, left = start, right = stop)))
}

map_df_reg <- function(pop, df){
    dplyr::select(df, -region_id) %>%
        purrr::pmap(read_df_reg, pop = pop)
}


merge_df <- function(snp_df, gwas_df){

  mutate(gwas_df, match_id = ldmap::find_alleles(chrom,
                                                 pos,
                                                 ref_chrom = snp_df$chr,
                                                 ref_pos = snp_df$pos)) %>%
      filter(!is.na(match_id)) %>%
      mutate(
          flip_allele = ldmap:::flip_alleles(
                                    allele,
                                    target_ref_alt = snp_df$allele[match_id]),
          snp_id = snp_df$snp_id[match_id]) %>%
      filter(flip_allele != 0) %>%
      mutate(`z-stat` = `z-stat` * flip_allele) %>%
      dplyr::select(-match_id, -flip_allele)
}

map_merge_df <- function(snp_df_l, gwas_df_l){
    purrr::map2(snp_df_l, gwas_df_l, merge_df)
}

merge_map <- function(inp_df, map_df){

    if (nrow(map_df) > 2){
        ret <- mutate(inp_df,
                      map = ldmap::interpolate_genetic_map(map = map_df$map,
                                                           map_pos = map_df$pos,
                                                           target_pos = pos,
                                                           strict = F))
        return(ret)
    }
    return(NULL)
}

map_merge_map <- function(inp_df_l, map_df_l){
    purrr::map2(inp_df_l, map_df_l, merge_map)
}


local_quh_gen <- function(inp_df, pop){
    if (is.null(inp_df)){
        return(NULL)
    }
    chrom <- unique(inp_df$chrom)
    stopifnot(length(chrom) == 1)
    input_file <- fs::path(kg_dir, glue("{pop}.chr{chrom}.h5"))
    X <- t(EigenH5::read_matrix_h5v(input_file, "dosage", inp_df$snp_id))
    gc()
    evd_R <- eigen(
        ldshrink::ldshrink(
                      genotype_panel = X,
                      map_data = inp_df$map,
                      na.rm = FALSE))
    gc()

    inp_df %>%
        dplyr::mutate(quh = c(t(evd_R$vectors) %*% `z-stat`),
                      D = evd_R$values)

}

map_local_quh_gen <- function(inp_df_l, pop){
    purrr::map(inp_df_l,
               local_quh_gen,
               pop = pop)
}

local_rssp_est <- function(o_df){

    out <- RSSp::RSSp_estimate(quh = o_df$quh,
                               D = o_df$D,
                               sample_size = mean(o_df$N))
    dplyr::mutate(out,
                  chr = o_df$chrom[1],
                  start = min(o_df$pos),
                  stop = max(o_df$pos))

}

map_local_rssp_est <- function(o_df_l){
    purrr::map(o_df_l, local_rssp_est)
}

rba <- function(...){
    list(...) %>%
        purrr::map_df(flatten_dfr)
}
