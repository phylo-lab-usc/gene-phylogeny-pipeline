## Introduction

### *By Fiel Dimayacyac and Doris Wu* 

### Background
Biologists pioneering the novel field of comparative functional genomics (CFG) attempt to infer the mechanisms of gene regulation by looking for similarities and differences of gene activity over time across multiple species [[1](#references)]. CFG research has provided important insights into functional attributes of the genome such as expression, chromatin accessibility, and transcriptional regulation. Through comparative studies, we are now able to investigate functional genomics across species and its role in the differentiated processes of development, morphology, physiology, and other phenotypes[[2](#references)]. 

Usually, three types of data are used as input for the analysis: functional data such as gene activity measurements, pathway data that represent a series of reactions within a cellular process, and phylogenetic relationship data that describe the relatedness of species[[3](#references)]. However, this field has faced various statistical challenges due to the absence of standardized tools or models used for visualizing all aspects of this comparative functional genomics dataset. The current CFG field has also failed to take phylogenetic approaches and evolutionary relationships across species into consideration[[4](#references)][[5](#references)]. 

To address the aforementioned challenges and opportunities, the Pennell Lab works towards developing new phylogenetic comparative functional genomics (pCFG) approaches to establish a firm statistical and analytical foundation for later CFG studies. In the current project, our research aim is to assess and develop models of evolution and phylogenetic comparative mothods for CFG. To do this, we use Arbutus to assess three different evolutionary models - Brownian Motion (BM),Ornstein-Uhlenbeck (OU) and EB - in terms of the model fit and model adequacy. Arbutus is a package developed to assess the adequacy of continuous trait models, more specifically, it assesses the fit of evolutionary models in an absolute sense[[6](#references)]. The analysis results would suggest how adequate the models currently employed in the field in terms of answewing the questions we are trying to ask. 


We also extend our preliminary analyses to compare gene family trees and species trees in terms of model fit and adequacy. According to our hypothesis, these two types of phylogenies are likely different due to the different evolutionary relatioships they portray. Species phylogenies are built based upon the overall genetic characteristics, therefore, despite its ability to describe the overall relationship, it fails to accurately depict each specific gene or gene family's evolutionary relationship [[7](#references)][[8](#references)]. In contrast, gene family phylogenies are more accurate as they represent the evolutionary history of the gene studied. Therefore, here, we aim to investigate this hypothesis by assessing the model fit and adequacy when gene family phylogenies are used (in comparison to when species phylogenies are used).


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
![image](https://user-images.githubusercontent.com/92889727/163547073-83c7a7b5-b1a8-4c73-81de-b4a25d37c442.png)


1. Data Acquisition: The data used in the current analysis was obtained from Dr. Cope in his recent study[[9](#references)]. 
   
2. Converting script: a R script designed to convert the CDS files with protein sequences grouped by species into files grouped by gene families. This step is necessary for later alignment and tree building for gene familie phylogenies. In this script, getGene functions are specific to each species for gene ID isolation. Subsequently, all the data will be loaded into a single dataframe. The switch function uses a specific function depending on a condition. In this case, we are using the name of the file/species. Subsequently, the script uses single data frame with all the data in inner joins to separate by gene ortholog. The orthologs matrix and inner join is flipped so that each row represents a gene famiy. The gene sequences are then left joined by each row.

3. Alignment using MAFFT[[10](#references)]: This is the first step of the snakemake pipeline. The CDS fasta files of each gene family is used as input and the sequences are aligned together using MAFFT (a program for multiple sequence alignment). 

4. FastTree[[11](#references)]: From the aligned sequences, approximately-maximum-likelihood phylogenetic trees are built using FastTree.   
5. Trimming: As shown in the ortholog matrix, for a speicfic gene family, some species may have missing gene expression data. In this case, the species without expression data will be trimmed off the gene phylogenetic tree.

5. Treeparser: This step turns the trees generated into tree objects that are readble by R. This is necessary for using the trees as input for the subsequent arbutus analysis.

6. Arbutus analysis: This is the final step of the snakemake pipeline and it analyses and generates the results for investigation on the model fit and adequacy. 


### Dataset used 
Three files are used as inputs: *CDS.tar.xz* (containing .fasta files for genes in each species), *coevolution_data.csv* and *ortholog_final_matrix_0927.tsv* provided by Dr. Cope (cite). These can be found under the data directory. 

The CDS fasta files are organized according to species. Within each file(i.e. species), the genes are listed individually with a header describing the gene identity. The ortholog final matrix summarizes the ortholog IDs of the same genes in different species. The coevolution data matrix contains expression data of all the orthologs.

## Outputs and Result Analysis
The outputs should be in the /gene_phylogenies/ directory. Specifically, gene trees should be found in /gene_phylogenies/trees. With our input files, a total of 2319 gene family phylogenetic trees was formed in this process.

Using gene expression data and gene phylogenies, the arbutus analysis generates model fit and model adequacy results. 
First, the relative fit AIC results are generated: the fit of three different models (BM, OU and EB) is compared against each other to identify the model with the best fit. 
![image](https://user-images.githubusercontent.com/92889727/163086781-89cbdc28-f7e9-4dd1-b343-978caa515326.png)

As shown in the AIC figure above, Ornstein-Uhlenbeck (OU) is the model with the best fit for the majority of the genes. 

Then, to assess the adequacy of the best model (OU), six test statistics were chosen to balance statistical intuition and computational effort: 
1. *c.var*: the coefficient of variation (standard deviation/mean) of the absolute value of the contrasts
2. *d.cdf*: the D statistic obtained from a Kolmolgorov-Smirnov test from comparing the distribution of contrasts to that of a normal distribution with mean 0 and standard deviation equal to the root of the mean of squared contrasts
3. *m.sig*: the mean of the squared contrasts
4. *s.asr*: the slope of a linear model fitted to the absolute value of the contrasts against the ancestral state inferred at the corresponding node
5. *s.hgt*: the slope of a linear model fitted to the absolute value of the contrasts against node depth
6. *s.var*: the slope of a linear model fitted to the absolute value of the contrasts against their expected variances
![image](https://user-images.githubusercontent.com/92889727/163417655-d244ea4d-9eda-42a0-8e8c-6d8636bc057f.png)
All of these test statistics essentially evaluate whether the contrasts come from the distribution expected under BM. Simply speaking, we consider the model to be accurate if the test statistics are well distributed (i.e. nearly rectangular shaped).

----------------------

## References

[1] Lappalainen, T. (2015). Functional genomics bridges the gap between quantitative genetics and molecular biology. Genome Research, 35(10): 1427-1431. https://genome.cshlp.org/content/25/10/1427.full

[2] Brawand, D., Soumillon, M., Necsulea, A. et al. (2011). The evolution of gene expression levels in mammalian organs. Nature, 478: 343–348. https://doi.org/10.1038/nature10532

[3] Meyer, M., Wong, B., Munzner, T. et al. (2010). Pathline: A Tool For Comparative Functional Genomics. Computer Graphics Forum, 29(3): 1043-1052.  https://doi.org/10.1111/j.1467-8659.2009.01710.x

[4] Blake LE, Roux J, Hernando-Herraez I, Banovich NE, Perez RG, Hsiao CJ, Eres I, Chavarria C, Marques-Bonet T, Gilad Y (2018) A comparison of gene expression and DNA methylation patterns across tissues and species. bioRxiv, :487413. https://doi.org/10.1101/487413

[5] Cain CE, Blekhman R, Marioni JC, Gilad Y (2011) Gene expression differences among primates are associated with changes in a histone epigenetic modification. Genetics, 187(4):1225–1234.

[6] Pennell, M.W., Fitzjohn, R.G., Cornwell, W.K., & Harmon, L.J. (2015). Model Adequacy and the Macroevolution of Angiosperm Functional Traits. The American Naturalist, 186(2) https://doi.org/10.1086/682022

[7] Maddison, W., & Knowles, L.P. (2006). Inferring Phylogeny Despite Incomplete Lineage Sorting. Systematic Biology, 55(1): 21-30. https://doi.org/10.1080/10635150500354928

[8] Liu, L., & Pearl, D.K. (2007). Species Trees from Gene Trees: Reconstructing Bayesian Posterior Distributions of a Species Phylogeny Using Estimated Gene Tree Distributions. Systematic Biology, 56(3): 504-514. https://doi.org/10.1080/10635150701429982
 
[9] Cope, A.L., O'Meara, B.C., & Gilchrist, M.A. (2020). Gene expression of functionally-related genes coevolves across fungal species: detecting coevolution of gene expression using phylogenetic comparative methods. BMC Genomics, 21(370). https://doi.org/10.1186/s12864-020-6761-3

[10] Rozewicki, J., Li , S., Amada, K.M., et al. (2019). MAFFT-DASH: integrated protein sequence and structural alignment. Nucleic Acids Research 47: W5-W10. https://doi.org/10.1093/nar/gkz342

[11] Price, M.N., Dehal, P.S., & Arkin, A.P. (2009). FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix. Molecular Biology and Evolution, 26(7): 1641-1650. https://doi.org/10.1093/molbev/msp077
