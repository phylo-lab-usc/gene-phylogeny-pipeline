#Snakefile

IDS, = glob_wildcards("gene_phylogenies/fasta_files/{id}.fa")

rule all:
  input: expand("gene_phylogenies/trees/{id}.tree", id = IDS)
  
rule align:
  input:
    "gene_phylogenies/fasta_files/{id}.fa"
  output:
    "gene_phylogenies/alignments/{id}"
  shell:
    "mafft --auto --thread 8 {input} > {output}"
    
rule trim:
  input:
    "gene_phylogenies/alignments/{id}"
  output:
    "gene_phylogenies/trimmed/{id}"
  shell:
    "trimal -in {input} -out {output}"
    
rule tree:
  input:
    "gene_phylogenies/trimmed/{id}"
  output:
    "gene_phylogenies/trees/{id}.tree"
  shell:
    "FastTree -nt {input} > {output}"