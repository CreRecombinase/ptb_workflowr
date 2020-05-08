read_susie <- function(model,L=1,p_rank=1,region_id,...){
  target_name <- glue::glue("susie_res_{p_rank}L_{L}L_prior_r_mix_results_.scratch.midway2.nwknoblauch.ptb_scratch.{model}.txt.gz")



  set_r <- pluck(target_res,"susie_res","sets","cs",.default = list())
  sets <- flatten_int(set_r)
  stopifnot(all(as.integer(target_res$df$region_id)==as.integer(region_id)))
  ret <- mutate(target_res$df,
         susie_id=1:dplyr::n(),
         credible_set=susie_id %in% sets,best_pip=pip==max(pip)) %>% 
    dplyr::select(-susie_id) %>% dplyr::left_join(null_model_df) %>% 
    mutate(seqnames=paste0("chr",chrom),log10_p=-log10(p),width=1) %>%
    dplyr::select(seqnames,start=pos,log10_p,prior,pip,credible_set,best_pip,null_pip,t_r,width,rsid=id) %>% 
    plyranges::as_granges() %>%
    plyranges::set_genome_info("hg19")

}
