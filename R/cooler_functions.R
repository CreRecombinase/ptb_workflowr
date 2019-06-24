bed2cool <- function(cp){

  bed_f <- fs::dir_ls("/home/nwknoblauch/Dropbox/scratch/ptb_scratch/new_bed/",glob = "*bed.bz2")



  atac_files_p <-bed_f[str_detect(bed_f,"atac-seq-pooled")]
  af <- bzfile(atac_files_p[1])
  tdf <- read_tsv(af,col_names=c("chrom","start","stop")) %>% mutate(value=runif(n())) %>%
  ifd <- rtracklayer::import(af,format="BED")
  of <- fs::path("../../scratch/ptb_scratch/coolers/",fs::path_ext_set(fs::path_ext_remove(fs::path_file(atac_files_p)),rep(".bedGraph",length(atac_files_p))))
  rtracklayer::export.bedGraph(ifd,of[1])
  id <- str_replace(atac_files_p,".+pooled-(DSC[123])-.+","\\1")
  trt <- str_replace(atac_files_p,".+pooled-DSC[123]-([a-z]+)-ATAC.+","\\1")

  atac_files_np<-bed_f[str_detect(bed_f,"atac-seq-pooled",negate = T)&str_detect(bed_f,"atac-seq")]




  cs <- "/home/nwknoblauch/Dropbox/scratch/ptb_scratch/hic/hg19.chrom.sizes"
  ocs <- "/home/nwknoblauch/Dropbox/scratch/ptb_scratch/hic/c_hg19.chrom.sizes"
  read_tsv(cs,col_names = c("chrom1","size")) %>%
    filter(chrom1 %in% glue::glue("chr{c(1:23,'X')}")) %>%
    write_tsv(ocs,col_names = F)
  tempfgz <- fs::file_temp(ext = ".tsv.gz")
  tcsv <- readr::read_tsv(cp) %>%
    select(chrom1=bait_chr,pos1=bait_start,chrom2=otherEnd_chr,pos2=otherEnd_start,data=N_reads,score) %>%
      filter(chrom1 %in% glue::glue("chr{c(1:23,'X')}"),
             chrom2 %in% glue::glue("chr{c(1:23,'X')}")) %>%
    write_tsv(tempfgz,col_names = F)

}

res <- httr::GET("http://localhost:8989/api/v1/tilesets/") %>% content("parsed") %>% roomba(cols="uuid")


nres <- httr::GET("http://localhost:8888/api/v1/tiles/",query=list(d=paste0(res$uuid[1],".0.0.0"))) %>% content("parsed")

