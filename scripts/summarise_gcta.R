library(purrr)
library(readr)
library(dplyr)

inf <- snakemake@input[["in_f"]]
fgeneid <- snakemake@params[["gene"]]
outf <- snakemake@output[["out_f"]]

stopifnot(!is.null(outf),!is.null(fgeneid),length(inf)==length(fgeneid))
map2_dfr(inf,fgeneid,~read_delim(.x,delim="\t") %>% mutate(fgeneid=.y)) %>% write_delim(outf,delim="\t")
