library(vroom)
library(purrr)
library(fs)

anno_bed_if <- snakemake@input[["anno_bed"]]
l2m_if <- snakemake@input[["l2"]]
l2m50_if <- snakemake@input[["l2M_50"]]


l2m_of <- snakemake@output[["l2"]]
l2m50_of <- snakemake@output[["l2M_50"]]
pull_features <- snakemake@params[["features"]]
anno_cols <- scan(anno_bed_if,what = character(),nlines = 1)[-c(1:4)]
keep_cols <- anno_cols %in% pull_features


l2md <- scan(l2m_if,what=character(),nlines = 1)[keep_cols]
l2m50d <- scan(l2m50_if,what=character(),nlines = 1)[keep_cols]

write(paste0(l2md,collapse = "\t"), l2m_of)
write(paste0(l2m50d,collapse = "\t"),l2m50_of)
