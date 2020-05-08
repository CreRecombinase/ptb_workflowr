library(ldmap)
library(dplyr)
library(purrr)

input_f <- snakemake@input[["bimf"]]
output_f <- snakemake@output[["snplistf"]]

walk2(input_f, output_f,
      function(input,output){
  read_plink_bim(input) %>%
    count(rsid) %>%
    filter(n==1) %>%
    pull(rsid) %>% 
    write_lines(output)
})
