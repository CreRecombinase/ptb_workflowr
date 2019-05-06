
snp_reg <- function(t_chrom,t_start,t_stop,beta_v,se_v,N_v){
  db_df <- input_db_f
  dplyr::tbl(dplyr::src_sqlite(path = db_df, create = F),"gwas") %>%
    dplyr::select(SNP = id,
                  chrom ,
                  pos,
                  MAJOR = A1,
                  MINOR = A2,
                  N= !!N_v,
                  beta = !!beta_v,
                  se = !!se_v)%>%
    dplyr::mutate(`z-stat` = as.numeric(beta) / as.numeric(se),N=as.numeric(N))%>%
    mutate(chrom = as.integer(chrom), pos = as.integer(pos))%>% dplyr::select(-beta,-se)%>%
    dplyr::filter(chrom == t_chrom,between(pos,t_start,t_stop)) %>%dplyr::collect() %>%
    tidyr::unite(col="allele",MAJOR,MINOR,sep=",") %>%
    dplyr::distinct(chrom,pos,.keep_all=T) %>%
    arrange(chrom, pos)
}


read_map <- function(chrom,start,stop){

  snp_df_c <- EigenH5::read_vector_h5(map_file,"SNPinfo/chr")

  which_c <-which(snp_df_c==chrom)
  stopifnot(length(which_c)>0)
  EigenH5::read_df_h5(map_file,"SNPinfo",subcols=c("pos","map"),subset=which_c) %>%
    filter(dplyr::between(pos,start,stop)) %>% distinct(map,.keep_all=T) %>% arrange(pos)
}

read_df_reg <- function(pop,chrom,start,stop,read_map=F){

  input_f <-fs::path(kg_dir,glue("{pop}.chr{chrom}.h5"))
  snp_df_p <- EigenH5::read_vector_h5(input_f,"SNPinfo/pos")
  subcols <-c("SNP","allele","chr","pos","snp_id")
  if(read_map){
  subcols <- c(subcols,"map")
  }
  EigenH5::read_df_h5(input_f,"SNPinfo",subcols=subcols,
             subset=which(dplyr::between(snp_df_p,left = start,right = stop)))
}

merge_df <-function(snp_df,gwas_df){

  mutate(gwas_df,match_id=ldshrink::find_alleles(chrom,
                                                 pos,ref_chrom =snp_df$chr,ref_pos = snp_df$pos)) %>%
    filter(!is.na(match_id)) %>%
    mutate(flip_allele=ldshrink:::flip_alleles(allele,target_ref_alt = snp_df$allele[match_id]),
           snp_id=snp_df$snp_id[match_id]) %>%
    filter(flip_allele!=0) %>% mutate(`z-stat`=`z-stat`*flip_allele) %>%
    select(-match_id,-flip_allele)
}

merge_map <- function(inp_df,map_df){

    if(nrow(map_df)>2){
        ret <- mutate(inp_df,map=ldshrink::interpolate_genetic_map(map = map_df$map,
                                                            map_pos = map_df$pos,
                                                            target_pos = pos,strict = F))
        return(ret)
    }
    return(NULL)
}


local_quh_gen <- function(inp_df, pop){
    if (is.null(inp_df)){
        return(NULL)
    }
    chrom <- unique(inp_df$chrom)
    stopifnot(length(chrom)==1)
    input_file <-fs::path(kg_dir,glue("{pop}.chr{chrom}.h5"))
    X <-t(EigenH5::read_matrix_h5v(input_file,"dosage",inp_df$snp_id))
    evd_R <- ldshrink::ldshrink_evd(reference_panel = X,map = inp_df$map)
    inp_df %>% dplyr::mutate(quh=c(t(evd_R$Q)%*%`z-stat`),D=evd_R$D)

}

local_rssp_est <- function(o_df){


  out <- RSSp::RSSp_estimate(quh = o_df$quh, D = o_df$D, sample_size = mean(o_df$N))
    dplyr::mutate(out, chr = o_df$chrom[1], start = min(o_df$pos), stop=max(o_df$pos))

}

