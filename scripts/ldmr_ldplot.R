library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ldmap)
library(EigenH5)
library(Matrix)

inputf <- snakemake@input[["ldhf"]]


saveRDS(filter(df,r2 > 0.1),snakemake@output[["dataf"]])
ld_plot  <- filter(df,r2 > 0.1) %>% ggplot() +
  geom_point(aes(positions(isnp) + y, y, color = r2, alpha = r2), size = rel(0.5)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "Position", y = NULL) +
  scale_alpha(guide = 'none')
saveRDS(ld_plot,snakemake@output[["plotf"]])
