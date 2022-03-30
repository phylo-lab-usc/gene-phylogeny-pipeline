#Load and separate data
library(tidyverse)
library(Biostrings)

#load final matrix
orthologs <- read.delim("data/ortholog_final_matrix_0927.tsv")

to_read <- list.files("data/CDS/")

#load one of the CDS datasets
afumigatus <- readDNAStringSet("data/CDS/afumigatus.fasta")
df <- data.frame(gene = names(afumigatus), seq = paste(afumigatus))

#Need a getGene function for every species
#This function chooses which of the split pieces is the gene ID. Use the orthologs object as a reference.
#For example, for afumigatus, genes are labeled as "AFUB_*****" according to the orthologs object (the file from line 6 above).
#I noticed that this was always the 10th entry in the character vector produced by using str_split, so I returned the 10th index in my function.
getGeneAfum <- function( fasta_df ){
  chr_vec <- fasta_df$gene %>% str_split("[: ]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[10]) %>% unlist()
  fasta_df
}

#The map function below just applies the getGeneAfum function to each object in the vector f
#After running this, check out the df object. Now it should nicely turn out to be gene names in one column and sequences in another.
df <- getGeneAfum(df)

#This is the getGene function for Anidulans. I noticed that the gene name was the first index in this case.
#Each species has a different naming convention, so make sure to make a getGene function for each species
getGeneAnid <- function( fasta_df ){
  chr_vec <- fasta_df$gene %>% str_split("[: ]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}

#This is the getGene function for C albicans
#This file has the gene ID after "locus_tag=CAALFM_"
getGeneCalb <- function( fasta_df ){
  chr_vec <- str_split(fasta_df$gene, "locus_tag=CAALFM_")
  fasta_df$gene <- lapply(chr_vec, function(x) x[2]) %>% unlist()
  chr_vec <- str_split(fasta_df$gene, "]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}

#This is the getGene function for C glabrata 
#This file has the gene ID after "locus_tag="
getGeneCgla <- function( fasta_df ){
  chr_vec <- str_split(fasta_df$gene, "locus_tag=")
  fasta_df$gene <- lapply(chr_vec, function(x) x[2]) %>% unlist()
  chr_vec <- str_split(fasta_df$gene, "]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}


#This is the getGene function for C neoformans 
#This file has the gene ID after "locus_tag="
getGeneCneo <- function( fasta_df ){
  chr_vec <- str_split(fasta_df$gene, "locus_tag=")
  fasta_df$gene <- lapply(chr_vec, function(x) x[2]) %>% unlist()
  chr_vec <- str_split(fasta_df$gene, "]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}

#This is the getGene function for C parapsilosis 
#This file has the gene ID after "locus_tag="
getGeneCpar <- function( fasta_df ){
  chr_vec <- str_split(fasta_df$gene, "locus_tag=")
  fasta_df$gene <- lapply(chr_vec, function(x) x[2]) %>% unlist()
  chr_vec <- str_split(fasta_df$gene, "]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}

#This is the getGene function for F graminearum 
#This file has the gene ID in the 10th field 
getGeneFgra <- function( fasta_df ){
  chr_vec <- fasta_df$gene %>% str_split("[: ]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[10]) %>% unlist()
  fasta_df
}

#This is the getGene function for Ikluyveri
#gene ID is the first field
getGeneLklu <- function( fasta_df ){
  chr_vec <- fasta_df$gene %>% str_split("[: ]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}

#This is the getGene function for M oryzae
# gene ID is after "locus_tag="
getGeneMory <- function( fasta_df ){
  chr_vec <- str_split(fasta_df$gene, "locus_tag=")
  fasta_df$gene <- lapply(chr_vec, function(x) x[2]) %>% unlist()
  chr_vec <- str_split(fasta_df$gene, "]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}

#This is the getGene function for n castelli
# gene ID is after "protein_id="
getGeneNcas <- function( fasta_df ){
  chr_vec <- str_split(fasta_df$gene, "protein_id=")
  fasta_df$gene <- lapply(chr_vec, function(x) x[2]) %>% unlist()
  chr_vec <- str_split(fasta_df$gene, "]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}

#This is the getGene function for N crassa
# gene ID is after "locus_tag="
getGeneNcra <- function( fasta_df ){
  chr_vec <- str_split(fasta_df$gene, "locus_tag=")
  fasta_df$gene <- lapply(chr_vec, function(x) x[2]) %>% unlist()
  chr_vec <- str_split(fasta_df$gene, "]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}


#This is the getGene function for N discreta
# gene ID is after "neudi1|" and after pulling out the number, we add "NEUDI_" to it
getGeneNdis <- function( fasta_df ){
  chr_vec <- fasta_df$gene %>% str_split("[|]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[3]) %>% unlist()
  fasta_df$gene <- paste("NEUDI",fasta_df$gene, sep="_")
  fasta_df
}

#This is the getGene function for N tetrasperma
# gene ID is after "locus_tag="
getGeneNtet <- function( fasta_df ){
  chr_vec <- str_split(fasta_df$gene, "locus_tag=")
  fasta_df$gene <- lapply(chr_vec, function(x) x[2]) %>% unlist()
  chr_vec <- str_split(fasta_df$gene, "]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}

#This is the getGene function for S bayanus
# gene ID is the first field
getGeneSbay <- function( fasta_df ){
  chr_vec <- fasta_df$gene %>% str_split("[: ]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}

#This is the getGene function for S cerevisiae
# gene ID is after "protein_ID="

getGeneScer <- function( fasta_df ){
  chr_vec <- str_split(fasta_df$gene, "protein_id=")
  fasta_df$gene <- lapply(chr_vec, function(x) x[2]) %>% unlist()
  chr_vec <- str_split(fasta_df$gene, "]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}

#This is the getGene function for S kudriavzevii
# gene ID is the first field
getGeneSkud <- function( fasta_df ){
  chr_vec <- fasta_df$gene %>% str_split("[: ]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}

#This is the getGene function for S mikatae
# gene ID is the first field
getGeneSmik <- function( fasta_df ){
  chr_vec <- fasta_df$gene %>% str_split("[: ]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}

#This is the getGene function for S paradoxus
# gene ID is the first field
getGeneSpar <- function( fasta_df ){
  chr_vec <- fasta_df$gene %>% str_split("[: ]")
  fasta_df$gene <- lapply(chr_vec, function(x) x[1]) %>% unlist()
  fasta_df
}


#This function isn't complete, but once all the getGene functions are made, we can use it to load all the data into a single data frame
#The switch function uses a specific function depending on a condition. In this case, we are using the name of the file/species!
read_and_fix <- function(name){
  path <- "data/CDS/"
  fas <- readDNAStringSet(paste0(path, name))
  df <- data.frame(gene = names(fas), seq = paste(fas))
  getGene <- switch(name,
                    afumigatus.fasta = getGeneAfum,
                    anidulans.fasta = getGeneAnid,
                    calbicans.fasta = getGeneCalb,
                    cglabrata.fasta = getGeneCgla,
                    cneoformans.fasta = getGeneCneo,
                    cparapsilosis.fasta = getGeneCpar,
                    fgraminearum.fasta = getGeneFgra,
                    lkluyveri.fasta = getGeneLklu,
                    moryzae.fasta = getGeneMory,
                    ncastellii.fasta = getGeneNcas,
                    ncrassa.fasta = getGeneNcra,
                    ndiscreta.fasta = getGeneNdis,
                    ntetrasperma.fasta = getGeneNtet,
                    sbayanus.fasta = getGeneSbay,
                    scerevisiae.fasta = getGeneScer,
                    skudriavzevii.fasta = getGeneSkud,
                    smikatae.fasta = getGeneSmik,
                    sparadoxus.fasta = getGeneSpar,
                    getGeneAfum)
  res <- getGene(df)
}

#Actual usage of read_and_fix function (incomplete)
big_data_frame <- map_df(to_read, read_and_fix)

#The next step would be to use the single data frame with all the data in inner joins to separate by gene ortholog
#Flip orthologs matrix and inner join by row. Because each row represents a gene family, I just need to left join by each row.
sep_into_gene_fams <- function(orthos){
  for(i in 1:nrow(orthos)){
    tmp1 <- orthos[i,] %>% t() %>% as.data.frame()
    tmp2 <- tmp1 %>%
      dplyr::rename(gene = as.character(i)) %>%
      left_join(big_data_frame, by = "gene") %>%
      drop_na()
    if(nrow(tmp2) > 1){
      name = tmp1[18,1]
      tmp3 <- tmp2 %>% select(seq, gene)
      saveRDS(tmp3, file = paste0("gene_phylogenies/tables/", name))
      biomaRt::exportFASTA(tmp3, file = paste0("gene_phylogenies/fasta_files/", name, ".fa"))
    }
  }
}

#An example of what each step in the function above does
#tmp <- orthologs[2,] %>% t() %>% as.data.frame() %>% dplyr::rename(gene = "2") %>% left_join(big_data_frame, by = "gene") %>%
#  drop_na() %>% select(seq, gene)
#biomaRt::exportFASTA(tmp, "test.fa")

#Running the function
sep_into_gene_fams(orthologs)