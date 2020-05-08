input_st
  inputf <- dir("/tmp",pattern="*bw",full.names=TRUE)
  susie_d <- "/home/nwknoblauch/Dropbox/Repos/ptb_workflowr/output/results/plot_data/susie"
  susie_f <- "/home/nwknoblauch/Dropbox/Repos/ptb_workflowr/output/results/plot_data/susie/ptb_treatedhic_512.RDS"
  susie_nullf <- "/home/nwknoblauch/Dropbox/Repos/ptb_workflowr/output/results/plot_data/susie/ptb_512.RDS"

  pip_df <- dplyr::inner_join(readRDS(susie_nullf),
                       readRDS(susie_f),by = c("snp_struct","pval"),suffix = c("_null","_model"))
  f_pip_df <- pip_df %>%  dplyr::arrange(snp_struct) %>% dplyr::mutate(lrmax = pmin(cummax(pmax(pip_model,pip_null)),rev(cummax(rev(pmax(pip_null,pip_model))))))
  fm_region <- dplyr::filter(f_pip_df,lrmax > 1e-6)  %>% dplyr::pull(snp_struct)  %>% convex_hull()
  region_r <- fm_region
  x <- inputf[1]

  library(rtracklayer)
  library(ldmap)
  library(fs)
  library(purrr)
  library(dplyr)
  library(stringr)
  library(plyranges)
  library(Gviz)


  library(biovizBase)
  library(Homo.sapiens)
  library(readr)


  #exp <- example_HiC()
  #inputf <- snakemake@input[["bigwigs"]]
  fnames <- fs::path_ext_remove((fs::path_file(inputf)))
  mark <-  str_remove(fnames,"DSC[0-9]-.+-")
  treatment <-  str_replace(fnames,"DSC[0-9]-(.+)-.+","\\1")
  (sample <-  str_replace(fnames,"(DSC[0-9])-(.+)-.+","\\1"))
  input_df <- tibble::tibble(file=inputf,name=fnames,mark=mark,treatment=treatment,sample=sample)
  ## region <- 356L
  ## region_r <- ldetect_EUR[region]
  library(biomaRt)
  bm <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

  ## dTrack2 <- Gviz::DataTrack(range=x, genome="hg19", type="l", chromosome=paste0("chr",chromosomes(region_r)), name="bigwig")


  ## Gviz::plotTracks(dTrack2,from=starts(region_r),to=ends(region_r))
  focus_gr <- function(ldmr){
    chrp  <- paste0("chr",chromosomes(ldmr))
    start <- starts(ldmr)
    end <- ends(ldmr)

    return(GenomicRanges::GRanges(seqnames = chrp,ranges=IRanges::IRanges(start = start,end = end)))
  }

  region_r <- ldetect_EUR[512]
  gr <- focus_gr(region_r)
  sn <- as.character(seqnames(gr))
  biomTrack <- BiomartGeneRegionTrack(genome="hg19", chromosome=sn,start=starts(region_r),end=ends(region_r),geneSymbols=TRUE,name="ENSEMBL", filter=list(with_refseq_mrna=TRUE),biomart=bm,stacking="squish")

  top2 <- filter(pip_df, pip_model == max(pip_model))


  igrf <- function(tdf,ty){
    tyd <- ty$treatment
    inpl <- purrr::pmap(tdf, function(file, name, sample, ...){
      rtracklayer::import(file, which = gr)
    }) %>% plyranges::bind_regions() %>% disjoin_regions(score=mean(score)) 
    inpl %>% tibble::as_tibble() %>%
      dplyr::rename({{tyd}} := score) %>% 
      plyranges::as_granges()
  }

  ngrf <- function(idf,y, ...){
    rrangel <- group_by(idf,treatment) %>%
      group_map(igrf)  %>%  purrr::reduce(join_overlap_inner)
    dtrack <- Gviz::DataTrack(rrangel,groups=c("dec","unt"),name=glue::glue("{y$mark}"),type="hist",window=-1,windowSize=1500)
  }

  idf <- dplyr::filter(input_df,mark==mark[1])
  tdf <- filter(idf,treatment==treatment[1])
  rangel <- group_by(input_df,mark) %>% group_map(ngrf)

  trackl <- c(biomTrack,rangel)
  strt <- positions(top2$snp_struct)
  grangel <- HighlightTrack(trackl,chromosome=sn,start=strt,width=10)
  png(file="/home/nwknoblauch/Dropbox/Hand2.png")
  plotTracks(c(grangel),chromosome=sn,from=starts(region_r),to=ends(region_r))
  dev.off()

  GeneRegionTrack(txdb)
  rangeld <- lapply(rangel,function(df){
    mrk <- df$mark[1]
    df <- plyranges::select(df,-mark)
 
  })

    rangeldf <- bind_regions(rangel)


  plotTracks(rangeld,groups=c("DSC1","DSC2","DSC3"),type=c("a","p"))






  as_grange_df <- function(df, ...) {
    rl <- rlang::list2(...)
    explode_ldmap_region(df) %>%
      dplyr::rename(seqnames = chrom) %>%
      mutate(seqnames = as.character(seqnames)) %>%
      mutate(!!!rl) %>%
      plyranges::as_granges()
  }

  as_ldmr <- function(df) {
    tibble::as_tibble(df) %>%
      dplyr::mutate(seqnames = as.integer(stringr::str_remove(seqnames,"chr"))) %>% 
dplyr::select(-width) %>% 
      compact_ldmap_region(chrom = "seqnames")
  }


  read_df <- as_ldmr(ranged) %>% dplyr::select(-strand)
  pread_df <- tidyr::pivot_wider(read_df,names_from=c("mark","sample","treatment"),values_from=c("score"))
  ldmdf <- as_ldmr(out.gr)
