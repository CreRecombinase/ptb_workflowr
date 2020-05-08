## save.image("/tmp/ig.RData")
## stop("whoa!")
## load("/tmp/ig.RData")
  library(dplyr)
  library(purrr)
  library(readr)
  library(vroom)
  library(EigenH5)
  library(ldmap)

  input_f <- snakemake@input[["inputf"]]
  index_f <-  snakemake@input[["indexf"]]
  chrom <- snakemake@params[["chrom"]]
  stopifnot(!is.null(chrom))
  schrom <- as.integer(chrom)
  output_f <- snakemake@output[["outputf"]]


  ind_spec <- cols_only(
    CHR = col_integer(),
    BP = col_double(),
    SNP = col_character()
  )

  gwas_type <- if_else(
    is.null(snakemake@params[["gwas_t"]]),
    "",
    paste0(".", snakemake@params[["gwas_t"]])
  )

beta_col <- glue::glue("beta{gwas_type}")
se_col <- glue::glue("se{gwas_type}")
N_col <- glue::glue("N{gwas_type}")
P_col <- glue::glue("pval{gwas_type}")

sel_cols <- c("snp_struct",
              beta_col,
              "A1",
              "A2",
                se_col,
                N_col,
                P_col)

  sel_cols <- stringr::str_replace(
                         sel_cols,
                         "\\.$",
                         "")

  index_df <- vroom::vroom(
                       index_f,
                       delim = "\t",
                       col_types = ind_spec
                     )  %>%
    rename(chrom = CHR, rsid = SNP, pos = BP) %>%
    mutate(snp_pos=new_ldmap_snp(chrom,pos)) %>%
    select(-chrom,-pos)

  nr_index_df <- nrow(index_df)


input_i <- read_snp_region_h5(input_f,
                              ldmr = hg19_sizes[schrom],
                              datapath = "snp",
                              subcols = sel_cols) %>%
  rename(gwas_snp = snp_struct) %>%
  mutate(snp_pos = clear_alleles(gwas_snp)) %>%
  inner_join(index_df)

stopifnot(all(input_i$chrom == schrom))
stopifnot(nrow(input_i)>0)

input_i  %>% rename(beta =  {{beta_col}},
                    se =  {{se_col}},
                    N =  {{N_col}}) %>%
  dplyr::distinct(rsid, .keep_all = TRUE) %>%
  dplyr::transmute(SNP = rsid,
                   N = N,
                   Z = beta / se,
                   A1 = A1,
                   A2 = A2,
                   P=pval) %>%
  vroom::vroom_write(output_f, delim = "\t")
