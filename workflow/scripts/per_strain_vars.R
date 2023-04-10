library(rtracklayer)
library(tidyverse)
library(VariantAnnotation)
library(BSgenome.Dmelanogaster.UCSC.dm6)

# This script takes tidal TE ins/dep annotation.txt files and
# produces Granges objects containing all information necessary for conversion
# to vcf.

# for adding ref alleles
g <- BSgenome.Dmelanogaster.UCSC.dm6
seqlevelsStyle(g) <- "NCBI"

#ins_fl <- "results/polymorphic_insertions/dl/DGRP_flies/RAL-101_result/RAL-101_Inserts_Annotated.txt"
#deps_fl <- "results/polymorphic_insertions/dl/DGRP_flies/RAL-101_result/RAL-101_Depletion_Annotated_TEonly.txt"

ins_fl <- snakemake@params[["poly_ins_tsv"]]
deps_fl <- snakemake@params[["poly_deps_tsv"]]

ins <- read_tsv(ins_fl)
deps <- read_tsv(deps_fl)

sample_name <- ifelse(exists("snakemake"), snakemake@params$sample_name, "DGRP_101")

sample_name <- str_replace(sample_name,"RAL-","DGRP_")

# -------------------------------- process ins and deps

# see http://www.bio.brandeis.edu/laulab/Tidal_Fly/UserGuide_TIDAL_outputs.html

ins <- ins %>% 
  mutate(strand=ifelse(TE_coord_start > TE_coord_end,"-","+")) %>%
  dplyr::rename(chr=Chr,start=Chr_coord_5p,end=Chr_coord_3p,name=TE) %>%
  mutate(start = start) %>%
  mutate(chr=str_remove(chr,"chr"))%>%
  mutate(source = "Tidal") %>%
  mutate(ref = ".", alt="<INS>") %>%
  mutate(sample_name = sample_name) %>%
  mutate(refDepth=round(Norm_RefGen_Reads),altDepth=Reads_collapsed,totalDepth=Reads_collapsed + refDepth) %>%
  dplyr::select(chr,start,end,strand,name,sample_name,refDepth, altDepth, totalDepth,ref,alt) %>%
  GRanges()

deps <- deps %>% 
  filter(Chr_5p...2 == Chr_3p...24) %>%
  filter(Chr_5p...2 == Chr_mid) %>%
  dplyr::rename(chr=Chr_5p...2,start=Chr_coord_5p_end,end=Chr_coord_3p_start,name=repName) %>%
  mutate(start_old = start, end_old = end, start=ifelse(start_old > end_old,end,start),end=ifelse(start_old > end_old,start_old,end)) %>%
  mutate(start = start) %>%
  mutate(source = "Tidal") %>%
  mutate(chr=str_remove(chr,"chr"))%>%
  mutate(ref = ".", alt="<DEL>") %>%
  mutate(strand = "*") %>%
  mutate(sample_name = sample_name) %>%
  mutate(refDepth=round(RefGen_Avg),altDepth=Reads_collapsed,totalDepth=Reads_collapsed + refDepth) %>%
  dplyr::select(chr,start,end,strand,name,sample_name,refDepth, altDepth, totalDepth,ref,alt) %>%
  GRanges()

# note: there are some very large dels in the tidal data. need to consider what to do with these.
# if they're real don't want to exclude them or otherwise we could find a polymorphic TE
# that looks like it affects gene expression, but really its just a deletion that encompasses a gene.


gr <- c(ins,deps) %>% sort()
strand(gr) <- "*"
gr$ref <- Biostrings::getSeq(g, resize(gr,fix = "start",width = 1),as.character=T)

# this can be output as granges obj for
# downstream conversion to vcf or VEP tab-delimited input format
saveRDS(gr,snakemake@output$rds)