library(vroom)
library(purrr)
library(fs)

anno_bed_if <- snakemake@input[["anno_bed"]]
l2gz_if <- snakemake@input[["l2gz"]]

anno_bed_of <- snakemake@output[["anno_bed"]]
l2gz_of <- snakemake@output[["l2gz"]]
pull_features <- snakemake@params[["features"]]

annot_prefix <- c("CHR","BP","SNP","CM")
annot_cols <- c(annot_prefix, pull_features)
names(annot_cols) <- annot_cols
annot_cols <- map(annot_cols, ~col_guess())
annot_cols <- vroom::cols_only(!!!annot_cols)

l2_prefix <- c("CHR","SNP","BP")
l2_cols <- c(l2_prefix, paste0(pull_features, "L2"))
names(l2_cols) <- l2_cols
l2_cols <- map(l2_cols, ~col_guess())
l2_cols <- vroom::cols_only(!!!l2_cols)

vroom::vroom_write(vroom::vroom(l2gz_if, delim = "\t",col_types = l2_cols),l2gz_of,delim = "\t")
vroom::vroom_write(vroom::vroom(anno_bed_if,delim = "\t",col_types = annot_cols),anno_bed_of,delim = "\t")
