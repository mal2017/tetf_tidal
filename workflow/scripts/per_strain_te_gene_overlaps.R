library(tidyverse)
library(rtracklayer)
library(plyranges)

# ----- process ins to grl ------
#ins_fl <- "results/collected_insertions/insertions_by_strain.tsv.gz"
ins_fl <- snakemake@input[["ins"]]

ins <- read_tsv(ins_fl)

#genes_fl <- "../tetf_refs/results/combined-anno/combined.gtf.gz"
genes_fl <- snakemake@input[["genes"]]

genes <- rtracklayer::import(genes_fl)

genes <- genes[genes$type == "gene"]

genes$name <- genes$gene_id

maxgap <- 5000
maxgap  <- snakemake@params[["maxgap"]]

# import and remove simple repeats
ins_gr <- GRanges(ins)
ins_gr <- ins_gr[!str_detect(ins_gr$name,"\\([ACTG]+\\)n")]

# find overlapping genes and TEs for each strain
res0 <- join_overlap_inner(genes,ins_gr, maxgap = maxgap) %>%
  as_tibble() %>%
  mutate(gene_symbol.y=name.y) %>%
  dplyr::select(Strain=strain,name.x, gene_symbol.x=gene_symbol,
                name.y, gene_symbol.y) %>%
  distinct() %>%
  mutate(overlap=T)

# downstream analyses expect that name.x and name.y hold 
# symmetrical information
res <- res0 %>%
  dplyr::rename(gene_symbol.x=gene_symbol.y,gene_symbol.y=gene_symbol.x,
                name.x=name.y, name.y=name.x) %>%
  bind_rows(res0,.)


write_tsv(res,snakemake@output[["tsv"]])
