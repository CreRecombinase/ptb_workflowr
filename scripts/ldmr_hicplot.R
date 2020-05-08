#  #source("renv/activate.R")

  library(stringr)
  library(readr)
  library(EigenH5)
  library(ldmap)
  library(tidyr)
  library(vroom)
  library(purrr)
  library(dplyr)

  anno_file <- snakemake@input[["anno_file"]]
  region <- snakemake@params[["region"]]
  hic_file <- snakemake@input[["hic"]]
  null_file <- snakemake@input[["null_file"]]
  alt_file <- snakemake@input[["alt_file"]]

  pip_df <- inner_join(readRDS(null_file),
                       readRDS(alt_file),
                       by = c("snp_struct","pval"),
                       suffix = c("_null","_model")) %>%
    mutate(region_id = region)

  tfr <- pip_df %>%
    dplyr::arrange(snp_struct) %>%
    dplyr::mutate(lrmax = pmin(cummax(pip_model),rev(cummax(rev(pip_model))))) %>%
     dplyr::filter(lrmax > 1e-6)  %>% dplyr::pull(snp_struct)  %>% convex_hull()

    cs_df <-  dplyr::filter(pip_df,CS_model==TRUE)

  anno_df <- readRDS(anno_file) %>%
      mutate(region_id = region) %>%
      dplyr::filter(is_region_in_region(ldmap_region, tfr))

  hic_df <- readRDS(hic_file) %>%
      dplyr::mutate(region_id = y) %>%
    dplyr::filter(is_region_in_region(bait, tfr), is_region_in_region(target, tfr))

  bait_df <- dplyr::select(hic_df, ldmap_region = bait,
                           region_id) %>%
    mutate(anno = "DSC_Treated_HiC")

  target_df <- dplyr::select(hic_df, ldmap_region = target,
                             region_id) %>%
    mutate(anno = "DSC_Treated_HiC")

  plot_df <- bind_rows(anno_df,
                       bait_df,
                       target_df) %>%
    mutate(ldmap_region = as_ldmap_region(ldmap_region))


  midpoints <- function(x) {
      stopifnot(inherits(x, "ldmap_region"))
      return(starts(x) + round((ends(x) - starts(x)) / 2))
  }
  widths <- function(x) {
      stopifnot(inherits(x, "ldmap_region"))
      return(ends(x) - starts(x))
  }

  hlinel <- map(positions(cd_df$snp_struct), ~ geom_vline(xintercept = .x, col = "red", alpha = 0.2))

  add_plots <- function(proto, plot_list) {
      purrr::reduce(plot_list, `+`, .init = proto)
  }

  hicplotf <- function(reg) {
      fl <- focus_l[[reg]]

      hml <- hlinemapl[[reg]]
      df <- focus_hic_l[[reg]]

      tannogr <- focus_full_anno_l[[reg]] %>% dplyr::mutate(
          feature = str_remove(anno, "^.+dec-"),
          feature = str_remove(feature, "DSC_Treated_"),
          feat_fact = factor(feature)
      )

      df <- dplyr::mutate(df, feat_fact = factor("HiC", levels = levels(tannogr$feat_fact)))

      annoplot <- ggplot() +
          geom_tile(
              data = tannogr,
              aes(
                  x = midpoints(ldmap_region),
                  y = feat_fact,
                  width = widths(ldmap_region)
              ),
              height = .5
          ) +
          geom_curve(data = df, aes(
              x = midpoints(bait),
              y = feat_fact,
              xend = midpoints(target),
              yend = feat_fact
          )) +
          ylab("Functional Genomic Annotation") +
          xlab("Position") +
          xlim(c(starts(fl), ends(fl)))
      add_plots(annoplot, hml)
  }

  genel <- map(havegene_list, function(reg) {
      fl <- focus_l[[reg]]
      x <- focus_grtm[[reg]]
      y <- hlinemapl[[reg]]
      add_plots(x@ggplot + ylab("Genes") + xlim(c(starts(fl), ends(fl))), y)
  })
  fmplots <- list(gene = genel, anno = hicplots, ld = lddfl_vl) %>% pmap(function(gene, anno, ld) {
      (anno / gene / ld)
  })
  plan(multiprocess)

  od <- "/home/nwknoblauch/Dropbox/MOD paper 1/Fine_Mapping/"
  future_walk(seq_along(fmplots), function(i) {
      plt <- fmplots[[i]]
      name <- havegene_list[i]
      output_png <- fs::path(name, ext = "png")
      ggsave(filename = output_png, plot = plt, path = od)
  })


    library(EigenH5)
    cols <- rev(grDevices::colorRampPalette(c("#41AB5D","#BDBDBD"),
                                    space = "Lab")(200))
    ld_files <- fs::dir_ls("~/tmp/LD/",glob="*h5")
    ld_region <- str_replace(ld_files,".+/[0-9]+_([0-9]+).h5","\\1")
    ldf_df <- tibble(inputf=ld_files,region=ld_region) %>% dplyr::filter(region %in% names(focus_l))
    lddfl <- pmap(ldf_df,function(inputf,region){
      tfoc <- focus_l[[region]]
      tmpdf <- read_df_h5(inputf,"snp")
      tmpR <- read_matrix_h5(inputf,"R")
      tmpT <- as(tmpR, "dgTMatrix")
      upper <- (tmpT@i <= tmpT@j)
      df <- tibble(
        i = tmpT@i[upper], 
        j = tmpT@j[upper],
        r2 = tmpT@x[upper]^2
      ) %>% 
        mutate(isnp = tmpdf$snp_struct[i + 1],jsnp = tmpdf$snp_struct[j + 1]) %>% filter(positions(i) <=  positions(j)) %>% 
        mutate(y = (positions(jsnp) - positions(isnp)) / 2) %>% 
        dplyr::filter(is_snp_in_region(isnp,tfoc),is_snp_in_region(jsnp,tfoc)) %>%
        ggplot() +
        geom_point(aes(positions(isnp) + y, y, color = r2, alpha = r2), size = rel(0.5)) +
        scale_color_gradientn(colours=cols) +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        labs(x = "Position", y = NULL) +
        scale_alpha(guide = 'none')
    }) %>% set_names(ldf_df$region)

      focus_gr_list <- map(focus_l,function(ldmr) {
          chrp  <- paste0("chr",chromosomes(ldmr))
          start <- starts(ldmr)
          end <- ends(ldmr)

          return(GenomicRanges::GRanges(seqnames = chrp,ranges=IRanges::IRanges(start = start,end = end)))
        })
        focus_tgdl <- map(focus_gr_list,function(gr) {
          suppressMessages(tg_df <-
                             OrganismDbi::selectByRanges(
                                            x = Homo.sapiens,
                                            ranges = gr,
                                            columns = c("SYMBOL", "TXNAME")) %>%
                             as_tibble() %>%
                             dplyr::select(tx_name = TXNAME, symbol = SYMBOL) %>%
                             tidyr::unnest(cols = c(tx_name, symbol)))
          return(tg_df)
        })

      genefun <- function(tg_df,gr) {
        if(nrow(tg_df)>0){
          suppressMessages(gr.txdb <-
                             crunch(txdb, which = gr) %>%
                             tibble::as_tibble() %>%
                             dplyr::mutate(tx_name = as.character(tx_name)) %>%
                             dplyr::inner_join(tg_df) %>%
                             plyranges::as_granges() %>%
                             split(.$symbol))
        }else{
          gr.txdb <- GenomicRangesList()
          return(gr.txdb)
        }
      }


      focus_grtl <- map2(focus_tgdl, focus_gr_list, genefun)  %>% purrr::discard(~class(.x)=="SimpleGenomicRangesList")
      focus_grtm <- map(focus_grtl,autoplot)
      length(focus_grtl)

    txdb <- TxDb(Homo.sapiens)
