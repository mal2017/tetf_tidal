library(tidyverse)

#fl <- "results/overlaps/overlaps.tsv.gz"
fl <- snakemake@input$tsv

dat <- read_tsv(fl)

n_strains <- dat$Strain %>% unique() %>% length()

# get list of genes/te pairs that overlap
# in N strains out of the N DGRP strains in TIDAL
res <- dat %>%
  group_by(name.x, gene_symbol.x, name.y, gene_symbol.y) %>%
  summarise(n=sum(overlap),.groups = "drop") %>%
  filter(n==n_strains) %>%
  dplyr::select(-n)


write_tsv(res, snakemake@output$tsv)