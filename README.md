## Introduction

### *By Fiel Dimayacyac and Doris Wu* 

### Background
Biologists pioneering the novel field of comparative functional genomics (CFG) attempt to infer the mechanisms of gene regulation by looking for similarities and differences of gene activity over time across multiple species. CFG research has provided important insights into functional attributes of the genome such as expression, chromatin accessibility, and transcriptional regulation. Through comparative studies, we are now able to investigate functional genomics across species and its role in the differentiated processes of development, morphology, physiology, and other phenotypes. 

Usually, three types of data are used as input for the analysis: functional data such as gene activity measurements, pathway data that represent a series of reactions within a cellular process, and phylogenetic relationship data that describe the relatedness of species. However, this field has faced various statistical challenges due to the absence of standardized tools or models used for visualizing all aspects of this comparative functional genomics dataset. The current CFG field has also failed to take phylogenetic approaches and evolutionary relationships across species into consideration. 

To address the aforementioned challenges and opportunities, the Pennell Lab works towards developing new phylogenetic comparative functional genomics (pCFG) approaches to establish a firm statistical and analytical foundation for later CFG studies. In the current project, our research aim is to assess and develop models of evolution and phylogenetic comparative mothods for CFG. To do this, we use Arbutus to assess three different evolutionary models - Brownian Motion (BM),Ornstein-Uhlenbeck (OU) and EB - in terms of the model fit and model adequacy. We also extend our preliminary analyses to compare gene family trees and species trees in terms of model fit and adequacy. 


### Purpose
To achieve the aforementioned research aims, this pipeline was made to build gene phylognies given either protein accessions or sequence data. After the gene phylogenies are built, they can be used as input (along with gene expression data) for Arbutus. The model fit and adequacy analysis results will then be compared against the results for species phylogenies.

### Rationale 
By identifying the potential difference in model fit and adequacy when different phylogenies (specie versus gene families) are used, we will gain more insights into phylogenetic approaches and evolutionary models to use in CFG studies.


## Usage

### Installation

Installing this pipeline requires `conda` and `git`. Instructions for installing these
software can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) respectively. 

First, clone the repository by running the following command in a terminal:

```
git clone https://github.com/pennell-lab-ubc/gene-phylogeny-pipeline

```
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

1. Data Acquisition: The data used in the current analysis was obtained from Dr. Cope in his recent study (https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6761-3)
   
2. Converting script: a R script designed to convert the CDS files with protein sequences grouped by species into files grouped by gene families. This step is necessary for later alignment and tree building for gene familie phylogenies. In this script, getGene functions are specific to each species for gene ID isolation. Subsequently, all the data will be loaded into a single dataframe. The switch function uses a specific function depending on a condition. In this case, we are using the name of the file/species. Subsequently, the script uses single data frame with all the data in inner joins to separate by gene ortholog. The orthologs matrix and inner join is flipped so that each row represents a gene famiy. The gene sequences are then left joined by each row.

3. Alignment using MAFFT: This is the first step of the snakemake pipeline. The CDS fasta files of each gene family is used as input and the sequences are aligned together using MAFFT (a program for multiple sequence alignment). 

4. Trimming: As shown in the ortholog matrix, for a speicfic gene family, some species may have missing gene expression data. In this case, the species without expression data will be trimmed off the gene phylogenetic tree.

5. Treeparser: This step turns the trees generated into tree objects that are readble by R. This is necessary for using the trees as input for the subsequent arbutus analysis.

6. Arbutus analysis: This is the final step of the snakemake pipeline and it analyses and generates the results for investigation on the model fit and adequacy. 


### Dataset used 
Three files are used as inputs: *CDS.tar.xz* (containing .fasta files for genes in each species), *coevolution_data.csv* and *ortholog_final_matrix_0927.tsv* provided by Dr. Cope (cite). These can be found under the data directory. 

The CDS fasta files are organized according to species. Within each file(i.e. species), the genes are listed individually with a header describing the gene identity. The ortholog final matrix summarizes the ortholog IDs of the same genes in different species. The coevolution data matrix contains expression data of all the orthologs.

## Outputs and Result Analysis
The outputs should be in the /gene_phylogenies/ directory. Specifically, gene trees should be found in /gene_phylogenies/trees.  

Using gene expression data and gene phylogenies, the arbutus analysis generates model fit and model adequacy results. 
First, the relative fit AIC results are generated: the fit of three different models (BM, OU and EB) is compared against each other to identify the model with the best fit. Then, for that best model, we analyze the model fit and adequacy in terms of five aspects: ...

Second, we have arbutus results demonstrating the absolute fit. From the test statistics distribution, we can get an idea of how adequate and well distributed the results of a certain model are.


