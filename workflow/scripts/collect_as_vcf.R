# ---------------------------------------------------
library(VariantAnnotation)
library(tidyverse)
library(plyranges)
library(BSgenome.Dmelanogaster.UCSC.dm6)

# for adding ref alleles
g <- BSgenome.Dmelanogaster.UCSC.dm6
seqlevelsStyle(g) <- "NCBI"


fls <- Sys.glob("results/vcf/DGRP_flies/*rds")

gr <- lapply(fls,readRDS) %>%
  GRangesList() %>%
  unlist()

if (str_detect(snakemake@wildcards[["tidal_group"]],"DGRP")) {
  gr <- mutate(gr,strain  = str_extract(sample_name,"DGRP_\\d+"))
}

#gr$strain %>% unique()

# where multiple libraries exist for a strain, collapse each ins/dep to 1 entry per strain
gr <- gr %>% 
  sort() %>%
  group_by(name,alt,ref,strain) %>%
  reduce_ranges(refDepth = sum(refDepth), altDepth = sum(altDepth), totalDepth = sum(totalDepth))

# filter out really long svs. some of these make it into the TE dels annotated file from tidal,
# but on occasion are duplicated or have nested TE deletions within (though not all TEs contained by a large SV are
# always annotated as depleted, leading to some confusion...)
# here I remove the really large dels because we're interested in TEs
gr <- gr[(width(gr) < 15000),]


# make a standardized set of loci to avoid counting the following situation 
# as two separate insertions
# strainA -------------
# strainB  -------------
# these might have coordinates off by just a few bp, but for
# sake of simplicity we want to consider them as the same insertion
standardized <- gr %>%
  sort() %>%
  group_by(name,alt) %>%
  reduce_ranges() %>%
  mutate(.,ID = paste(name,str_remove_all(alt,"<|>"),str_replace_all(as.character(.),":|-","_"),sep="__"))


# join the standardized set (just used for ID)
# with the actual call set, retaining the ranges from the
# standardized set
gr.ID <- join_overlap_left(standardized, gr) %>%
  filter(name.x==name.y) %>%
  filter(alt.x==alt.y) %>%
  mutate(alt=alt.x,name=name.y) %>%
  dplyr::select(-alt.y, -name.y, -alt.x, -name.x)

gr.ID$ref <- Biostrings::getSeq(g, resize(gr.ID,fix = "start",width = 1),as.character=T)

stopifnot(length(gr.ID)==length(gr))

#vr_list <- gr.ID %>%
#  split(.,.$strain) %>%
#  map(~makeVRangesFromGRanges(.x,sampleNames.field = "strain")) 

vr <- makeVRangesFromGRanges(gr.ID,sampleNames.field = "strain")
#vr@softFilterMatrix <- FilterMatrix(matrix(rep(T,length(vr)),ncol = 1),filterRules = FilterRules(exprs = T))


writeVcf(vr,snakemake@params$vcf,index=T)
#writeVcf(vr_list$DGRP_101,"~/Downloads/DGRP_101.vcf",index=T)
#writeVcf(vr_list$DGRP_105,"~/Downloads/DGRP_105.vcf",index=T)

