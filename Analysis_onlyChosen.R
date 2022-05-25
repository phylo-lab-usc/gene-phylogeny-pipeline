#Arbutus Analysis

library(ape)
library(geiger)
library(phytools)
library(tidyverse)
library(arbutus)
library(parallel)

tree.file <- "gene_phylogenies/genetrees"
gene.exp <- "data/gene_expression_tpm_matrix_updated_Standard.LogNorm.tsv"

last4 <- function(string){
  substr(string,1,nchar(string)-5)
}

author_genes <- read_tsv(file = "data/control_set_consistent_w_bm.txt") %>%
  pivot_longer(cols = !c(Type, Score), values_to = "Prot") %>%
  select(!name) %>% unique() %>% pull(Prot)

tree_names <- list.files("gene_phylogenies/trees/") %>% as.list() %>% as.character() %>% map_chr(last4)
temp.df<- t(read.table(gene.exp,sep="\t",header=T,row.names = 1,stringsAsFactors = F))
exp.df <- temp.df[,author_genes]
col_names <- colnames(exp.df)[colnames(exp.df) %in% tree_names]
exp.df <- exp.df %>% as.data.frame() %>% dplyr::select(all_of(col_names)) %>% as.matrix()
trees <- readRDS(tree.file)


#replace gene IDs with species names
convert <- function(phylo){
  ndiscreta <- phylo$tip.label %>% str_which("NEUDI")
  cneoformans <- phylo$tip.label %>% str_which("CNAG")
  cparapsilosis <- phylo$tip.label %>% str_which("CPAR2")
  calbicans <- phylo$tip.label %>% str_which("C.......A")
  afumigauts <- phylo$tip.label %>% str_which("AFUB")
  spombe <- phylo$tip.label %>% str_which("SP.C")
  moryzae <- phylo$tip.label %>% str_which("MGG_")
  cglabrata <- phylo$tip.label %>% str_which("CAGL")
  ntetrasperma <- phylo$tip.label %>% str_which("NEUTE1DRAFT")
  fgraminearum <- phylo$tip.label %>% str_which("FGRAMPH1")
  anidulans <- phylo$tip.label %>% str_which("^AN")
  ncrassa <- phylo$tip.label %>% str_which("^NCU")
  lkluyveri <- phylo$tip.label %>% str_which("SAKL")
  skudriavzevii <- phylo$tip.label %>% str_which("Skud_")
  sbayanus <- phylo$tip.label %>% str_which("Sbay_")
  sparadoxus <- phylo$tip.label %>% str_which("Spar_")
  smikatae <- phylo$tip.label %>% str_which("Smik_")
  scerevisiae <- phylo$tip.label %>% str_which("NP_")
  ncastellii <- phylo$tip.label %>% str_which("XP_")
  phylo$tip.label[ndiscreta] <- "N.discreta"
  phylo$tip.label[cneoformans] <- "C.neoformans"
  phylo$tip.label[cparapsilosis] <- "C.parapsilosis"
  phylo$tip.label[calbicans] <- "C.albicans"
  phylo$tip.label[afumigauts] <- "A.fumigatus"
  phylo$tip.label[spombe] <- "S.pombe"
  phylo$tip.label[moryzae] <- "M.oryzae"
  phylo$tip.label[cglabrata] <- "C.glabrata"
  phylo$tip.label[ntetrasperma] <- "N.tetrasperma"
  phylo$tip.label[fgraminearum] <- "F.graminearum"
  phylo$tip.label[anidulans] <- "A.nidulans"
  phylo$tip.label[ncrassa] <- "N.crassa"
  phylo$tip.label[lkluyveri] <- "L.kluyveri"
  phylo$tip.label[skudriavzevii] <- "S.kudriavzevii"
  phylo$tip.label[sbayanus] <- "S.bayanus"
  phylo$tip.label[sparadoxus] <- "S.paradoxus"
  phylo$tip.label[smikatae] <- "S.mikatae"
  phylo$tip.label[scerevisiae] <- "S.cerevisiae"
  phylo$tip.label[ncastellii] <- "N.castellii"
  class(phylo) <- "phylo"
  phylo
}

#Fix gene trees
gene_names <- colnames(exp.df)

ret_prot <- function(phylo){
  res <- phylo$tip %>% str_subset("NP_........")
  res
}

temptrees <- trees %>% lapply(convert)

rename <- trees %>% lapply(ret_prot) %>% as.vector("character")

names(temptrees) <- rename
#Reorder and remove phylos that are null and have only 1 node
final_trees <- temptrees[gene_names]

rm(temptrees, rename, gene_names, col_names, trees)

#For each gene, trim tips off the phylogeny when that species is missing a gene
runFC <- function ( dat ){
  fitResults <- vector(mode = "list", length = ncol(dat))
  for(j in 1:ncol(dat)){
    df <- dat[,j] %>% as.data.frame() %>% drop_na() %>% as.matrix()
    tdf <- treedata(final_trees[[j]], df, sort = TRUE)
    phy <- tdf$phy
    data <- tdf$data
    fitBM <- fitContinuous(phy, data, model = "BM")
    fitOU <- fitContinuous(phy, data, model = "OU")
    fitEB <- fitContinuous(phy, data, model = "EB")
    aic <- c(fitBM$opt[["aic"]], fitOU$opt[["aic"]], fitEB$opt[["aic"]])
    fit <- ifelse(min(aic) == aic[1], list(c(fitBM, model = "BM")), 
                  ifelse(min(aic) == aic[2], list(c(fitOU, model = "OU")), 
                         list(c(fitEB, model = "EB"))))
    fitResults[j] <- fit
  }
  fitResults
}

model_count <- function (fit) {
  ou = 0
  bm = 0
  eb = 0
  for(f in fit){
    vec <- f
    ifelse(vec$model == "OU", ou <- ou + 1, ifelse(vec$model == "BM", bm <- bm + 1, eb <- eb + 1))
  }
  data.frame(OU = ou, BM = bm, EB = eb)
}

fitResults <- runFC(exp.df)
saveRDS(fitResults, "arbutus/fitResults_chosen")

df <- model_count(fitResults)

b <- df %>% pivot_longer(c(OU, BM, EB), names_to = "model")

b %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
ggsave("arbutus/AIC_chosen.png")

run_arb <- function (fits){
  class(fits) <- "gfit"
  arby <- tryCatch(arbutus(fits), error = function(f)NA)
}

arb_result <- mclapply(fitResults, run_arb, mc.cores = 16)
arb_result <- arb_result[!is.na(arb_result)]
pvals <- map_df(arb_result, pvalue_arbutus)
saveRDS(pvals, file = "arbutus/pvalues_df_chosen.rds")

p_piv <- pvals %>% pivot_longer(cols = everything(), names_to = "tstat")
saveRDS(p_piv, file = "arbutus/pvalues_table_chosen.rds")

p_piv %>% ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
ggsave("arbutus/arbutus_results_chosen.png")

p_inade <- pvals %>% select(!m.sig) %>% 
  transmute(c.less = c.var <= 0.05, sv.less = s.var <= 0.05, sa.less = s.asr <= 0.05, sh.less = s.hgt <= 0.05, d.less = d.cdf <= 0.05) %>% 
  transmute(inade = c.less + sv.less + sa.less + sh.less + d.less) %>% count(inade) %>% mutate(prop = n/sum(n)) %>% mutate(inade = as.character(inade))

p_inade %>% ggplot(aes(x = inade, y = n, fill = inade)) + geom_bar(stat = "identity") + 
  xlab("Number of inadequacies") + ylab("Number of genes") + ggtitle("Amount of genes by number of inadequacies") + theme_bw()
