#Parse Gene Trees into list

library(ape)
library(parallel)
library(tidyverse)
library(phytools)

tree_names <- list.files("gene_phylogenies/trees/") %>% as.list() %>% as.character()

chrono <- function( tree ){
  res <- chronos(tree, model = "clock")
}

parse_trees <- function ( treestring ){
  tree <- read.tree(file = paste0("gene_phylogenies/trees/", treestring))
  tree <- chrono(tree) %>% multi2di() %>% midpoint.root()
  tree
}

treelist <- mclapply(tree_names, parse_trees, mc.cores = 2) %>% discard(is.null) #%>% keep(~ Nnode(.x) > 1)

saveRDS(treelist, file = "gene_phylogenies/genetrees")
