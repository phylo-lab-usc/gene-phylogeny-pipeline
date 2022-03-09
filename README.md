# gene-phylogeny-pipeline
A pipeline to build gene phylogenies given either protein accessions or sequence data.

Note: Remember to activate the environment; `phyloenv` before using pipeline with `conda activate phyloenv` 

# OrthoFinder: phylogenetic orthology inference for comparative genomics

This pipeline incorporates OrthoFinder which is a fast, accurate and comprehensive platform for comparative genomics. OrthoFinder infers rooted gene trees for all orthogroups and identifies all of the gene duplication events in those gene trees.

1. To run its **help** text
You can run OrthoFinder: 'python OrthoFinder_source/orthofinder.py -h' or './OrthoFinder/orthofinder -h'. 
OrthoFinder should print its 'help' text.

2. To **run OrthoFinder on a directory of protein sequence fasta files**
'./OrthoFinder/orthofinder -f /OrthoFinder/ExampleData/'

2. To **run OrthoFinder on a directory of DNA sequence fasta files (in nucleotides)**
'./OrthoFinder/orthofinder -f /OrthoFinder/ExampleData/ -d'


For more detailed instructions on how to use OrthoFinder, see
Emms, D.M. and Kelly, S. (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biology 20:238

Emms, D.M. and Kelly, S. (2015) OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biology 16:157

