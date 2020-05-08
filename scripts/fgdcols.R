mc <- cols(
  chrom = col_factor(c(as.character(1:22), "X")),
  pos = col_double(), #Pos
  rsid = col_character(), #Rsid
  A1 = col_character(), #Effect_allele
  A2 = col_character(), #Non_effect_allele
  beta = col_double(), #Effect
  se = col_double(), #StdErr
  pval = col_double(), #P
  HetPVal = col_double(),
  N = col_double(),
  SNP = col_character()
)
data_delim <- " "
