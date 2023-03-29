library(tidyverse)
library(plyranges)
library(rtracklayer)

#fa_fl <- "../tetf_refs/results/flybase-anno/genome.fasta.gz"
fa_fl <- snakemake@input[["genome_fa"]]
fa <- import(fa_fl)

names(fa) <- names(fa) %>% str_extract(".+?(?=\\s)")

#insertions_path <- "results/collected_insertions/DGRP_flies.tsv"
insertions_path <- snakemake@input[["ins"]]

# import and ignore simple repeats
ins <- read_tsv(insertions_path)
ins <- GRanges(ins)
ins <- ins[!str_detect(ins$name,"\\([ACTG]+\\)n")]

ins <- ins %>% mutate(strain = str_extract(strain,"DGRP_\\d+"))

seqlengths(ins) <- seqlengths(fa)[seqlevels(ins)]

# Now get the reduced set and annotate with number
# of strains each insertion is found in,
# taking care to avoid double-qcounting overlaps from the same strain
# strain 1: ------       --------
# strain 2:      ----------
# should yield...
#           -------------------- penetrance = 2
# and not penetrance=3

tes_all <- ins %>%
  group_by(name,strain) %>%
  reduce_ranges() %>%
  ungroup()

tes_penetrance <- tes_all %>%
  group_by(name) %>%
  reduce_ranges(score = plyranges::n_distinct(strain)) %>%
  ungroup()

ins <- mutate(ins,name=paste(name,strain,source,sep=".")) %>% plyranges::select(name,score)

export(ins,snakemake@output[["all_ins"]],format="bb")
export(tes_penetrance,snakemake@output[["penetrance"]],format="bb")