# gene-phylogeny-pipeline

A pipeline to build gene phylogenies given either protein accessions or sequence data.

Before running the pipeline, make sure you have the correct environment. 

First, build the environment using the command below (in the terminal):

conda env create --file environment.yml

Next, ensure you activate the environment using the following terminal command:

conda activate genefams

Finally, to run the pipeline, make sure you are in the correct directory. Your working directory should be the same directory the Snakefile is in. 
If you are in the correct directory, run the command below:

snakemake -c6 

This will make snakemake use 6 cores. You can change the number of cores to however many you think is appropriate for this job. Make sure to keep the terminal running. You may want to do this in your terminal app ssh-ing into the server rather than through the Rstudio web. If you choose to do it in a terminal app instead, I would suggest doing it in a tmux. If you don't know what that is, let me know!

The outputs should be in the /gene_phylogenies/ directory. Specifically, gene trees should be found in /gene_phylogenies/trees.  