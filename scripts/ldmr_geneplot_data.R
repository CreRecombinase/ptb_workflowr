library(stringr)
library(Homo.sapiens)
library(RColorBrewer)

library(purrr)
library(dplyr)
library(biovizBase)
library(Homo.sapiens)
library(readr)

txdb <- TxDb(Homo.sapiens)

chrp <- snakemake@params[[1]][["chrom"]]
stopifnot(!is.null(chrp))

start <- snakemake@params[[1]][["start"]]
stopifnot(!is.null(start))
start <- as.integer(start)

end <- snakemake@params[[1]][["end"]]
stopifnot(!is.null(end))
end <- as.integer(end)

gr <- GenomicRanges::GRanges(seqnames = chrp,ranges=IRanges::IRanges(start = start,end = end))

suppressMessages(tg_df <-
                   OrganismDbi::selectByRanges(
                                  x = Homo.sapiens,
                                  ranges = gr,
                                  columns = c("SYMBOL", "TXNAME")) %>%
                   as_tibble() %>%
                   dplyr::select(tx_name = TXNAME, symbol = SYMBOL) %>%
                   tidyr::unnest(cols = c(tx_name, symbol)))

suppressMessages(gr.txdb <-
                   crunch(txdb, which = gr) %>%
                   as_tibble() %>%
                   mutate(tx_name = as.character(tx_name)) %>%
                   inner_join(tg_df) %>%
                   plyranges::as_granges() %>%
                   split(.$symbol))

dplyr::as_tibble(tg_df) %>% distinct(symbol) %>% write_tsv(snakemake@output[["genelistf"]])
saveRDS(unlist(gr.txdb),snakemake@output[["outputf"]])
