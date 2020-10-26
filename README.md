# plasmid-pipeline

# Introduction
this is a plasmid analysis pipeline..

# Requirements
ruamel.yaml  
click  
snakemake

# Install
Clone package:

> git clone https://github.com/Wanli-HE/plasmid-pipeline.git  
> cd plasmid-pipeline

all scripts are under the folder

# Usage
1. all databastes and dependent tools are downloaded and installed automatically with default version by script:  
> python ~/plasmid-pipeline/bin/run_main.py downloading 
     
before running this command, you need to install anaconda or miniconda. If you haven’t done it already you need to configure conda with the bioconda-channel and the conda-forge channel. This are sources for packages beyond the default one:  
> conda config --add channels defaults  
> conda config --add channels bioconda  
> conda config --add channels conda-forge

For version of the tools and databases, you can using  
> python ~/plasmid-pipeline/bin/run_main.py downloading -h  
to check, and also can change the peremters to choose the version that you need. and some tools are downloaded by conda YAML file, in ~/plasmid-pipeline/envs, you can also change the version from there. 





author: Wanli HE (wanli.he@bio.ku.dk)

cooperator: Franziska Klincke (franziska.klincke@bio.ku.dk)
            Joseph Nesme (joseph.nesme@bio.ku.dk)
            Søren Johannes Sørensen (sjs@bio.ku.dk)
