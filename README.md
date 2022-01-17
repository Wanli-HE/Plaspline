# Plaspline, a pipeline for plasmid analysis form metagenome

# Introduction
Nowadays, Numerous tools are developed based on different strategies for detecting plasmid from metagenome. Apparently, all of them have merits, but the demerits are also non-negligible. Thus, how to estimate them and how to combine the advantages of these tools are necessary. Besides, it is lack of a proven, systematic, and comprehensive workflow for theses multiple preprocessing and analytical steps, which are not only including detecting plasmid, but also including significant downstream analysis of plasmid.   
  
we are benchmarking these tools and strategies, building up *Plaspline* based on the benchmarked tools, which is trying to combine the advantages of tools for a better result of detecting plasmid. Moreover, multiple downstream analysis is also added into this workflow includes quality control, assembly, circular/linear plasmids isolation, circular/linear plasmids genome verification and classification, and plasmids-relative-genes analysis, which aims to get a comprehensive analysis of plasmid, both in gene and plasmid community level. 


# Requirements
ruamel.yaml=0.16.12  
click=7.1.2  
snakemake=5.25.0  
biopython=1.78  
python=3.9

# Install  
create a conda environment and insatll dependencies:  
```
> conda create -n plaspline -c conda-forge -c bioconda python=3.9 biopython=1.78 snakemake=5.25.0 click=7.1.2 ruamel.yaml=0.16.12
```

Clone package:  
```
> git clone https://github.com/Wanli-HE/Plaspline.git 
```  
all scripts are under the 'Plaspline' folder


# Downlanding dependent databases and tools:

*** make sure that you have enough space in your folder for it (50G maybe!!!!)

```
> python ~/Plaspline/bin/run_main.py downloading       
```
 

# Downloading example samples:

```
> wget 
```


# Usage




# Contact

## author:   
   Wanli HE (wanli.he@bio.ku.dk)

## cooperator:  
   Franziska Klincke (franziska.klincke@bio.ku.dk)  
   Joseph Nesme (joseph.nesme@bio.ku.dk)  
   Søren Johannes Sørensen (sjs@bio.ku.dk)
