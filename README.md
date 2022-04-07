## Introduction

### *By Fiel Dimayacyac and Doris Wu* 

### Background
Biologists pioneering the novel field of comparative functional genomics (CFG) attempt to infer the mechanisms of gene regulation by looking for similarities and differences of gene activity over time across multiple species. CFG research has provided important insights into functional attributes of the genome such as expression, chromatin accessibility, and transcriptional regulation. Through comparative studies, we are now able to investigate functional genomics across species and its role in the differentiated processes of development, morphology, physiology, and other phenotypes. Usually, three types of data are used as input for the analysis: functional data such as gene activity measurements, pathway data that represent a series of reactions within a cellular process, and phylogenetic relationship data that describe the relatedness of species. However, this field has faced various statistical challenges due to the absence of standardized tools or models used for visualizing all aspects of this comparative functional genomics dataset. The current CFG field has also failed to take phylogenetic approaches and evolutionary relationships across species into consideration. 

To address the aforementioned challenges and opportunities, the Pennell Lab works towards developing new phylogenetic comparative functional genomics (pCFG) approaches to establish a firm statistical and analytical foundation for later CFG studies. In the current project, our research aim is to assess and develop models of evolution and phylogenetic comparative mothods for CFG. To do this, we use Arbutus to assess three different evolutionary models - Brownian Motion (BM),Ornstein-Uhlenbeck (OU) and EB - in terms of the model fit and model adequacy. We also extend our preliminary analyses to compare gene family trees and species trees in terms of model fit and adequacy. 


### Purpose
To achieve the aforementioned research aims, this pipeline was made to build gene phylognies given either protein accessions or sequence data. After the gene phylogenies are built, they can be used as input (along with gene expression data) for Arbutus. The model fit and adequacy analysis results will then be compared against the results for species phylogenies.

### Rationale 
By identifying the potential difference in model fit and adequacy when different phylogenies (specie versus gene families) are used, we will gain more insights into phylogenetic approaches and evolutionary models to use in CFG studies.


## Usage

### Installation

Remember to build the environment using the command below in the terminal:
```
conda env create --file environment.yml
```

### Running the pipeline 

Before running the pipeline, ensure you activate the environment we created (see above) using the terminal command below:
```
conda activate genefams
```

Finally, to run the pipeline, make sure you are in the correct directory. Your working directory should be the same directory the Snakefile is in. 
If you are in the correct directory, run the command below:
```
snakemake -c6 
```
Note: the "-c6" here simply means this pipeline will be ran using six cores. You can change the number of cores to however many you think is appropriate for this job.

You can choose to have your terminal app ssh-ing into the server or just run the pipeline through the Rstudio server website.


### Pipeline overview
![image](https://user-images.githubusercontent.com/92889727/162136287-d11cab09-8563-453c-a581-3e7401fed4ab.png)

1. Data Acquisition
2. Converting script
3. Alignment using MAFFT
4. Trimming
5. Treeparser
6. Arbutus analysis 



### Dataset used 
Two files are used as inputs. The CDS files and the ortholog_final_matrix provided by Dr. Cope (cite).

## Outputs and Result Analysis
The outputs should be in the /gene_phylogenies/ directory. Specifically, gene trees should be found in /gene_phylogenies/trees.  

