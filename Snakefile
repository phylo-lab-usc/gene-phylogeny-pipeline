#Initializing

#rule all:
  #input:
    #"arbutus/arbutus.png"

#rule format:
  #shell:
    #"for f in *fa ; do python ~/xxx/OrthoFinder/tools/primary_transcript.py $f ; done"
    
#rule orthogroup:
  ##input: #Might not be necessary
    ##sequences = "data/CDS/{species}.fa"
    ##orthologs = "data/ortholog_final_matrix_0927.tsv"
  #output: 
    #"proteomes/xxx/Orthofinder/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv"
  #shell:
    #"orthofinder -f primary_transcripts/"
    
#rule ortho_convert:
  #input:
    #dictionary = "data/ortholog_final_matrix_0927.tsv"
    #orthogroups = "proteomes/xxx/Orthofinder/Comparative_Genomics_Statistics/N0.tsv"
  #output:
    #"arbutus_ready_matrix.tsv"
  #shell:
    #"some_script"

#rule run_R:
  #input:
    #trees = "xxx/GeneTrees"
    #expression_data = "data/coevolution_data.csv"
  #output:
    "arbutus/arbutus.png"
  #script:
    "scripts/analysis.R"