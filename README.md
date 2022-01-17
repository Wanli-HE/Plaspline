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
## prepare config.yaml and sample.json file.

```
> python  <~/Plaspline/bin/run_main.py> preprocessing --path <your reads folder> --adapter <your adapter file> --phix  <your phix file> --assembler <defult:'spades'> --skip_qc <defult:"False">
```
Notice: 
1. In your reads folder: plaspline can only recognize the raw reads (forward and resvers) file which suffix is "_R1" or "_r1" (only show forward in example).  
      ***must be gzip file. otherwise in "quality control step" it will raise error.

2. If you want to do the quality control step,  <your adapter file> and <your phix file> must be input;    if your input is qc reads, not necessary do the  quality control step,  you don`t need input <your adapter file> and <your phix file>, but the  "--skip_qc" must be "True".

3. the assembler defult is spades, you can also choose megahit as assembler, by input the params "--assembler megahit".

for more information, "python Plaspline/bin/run_main.py preprocessing -h" 
  
#Working
1. before true runing plaspline, it is better to check whether there is any command line errors, by dry run (" -n ").

```
  > python <~/Plaspline/bin/run_main.py> working  <"step"> -j 5  -n 
```
Notice:
  1. before running, make sure the cofig.yaml and samples.json are under the floder which you are running the command.
  
  2. the total threads that you are using is 5 (-j), in the "config.yaml" file also records the max threads of each step ("threads: 8" in config.yaml). so the total threads = (threads)*(job) = 5*8 =40. so make sure that you have enough threads in your system.
  

# Contact

## author:   
   Wanli HE (wanli.he@bio.ku.dk)

## cooperator:  
   Franziska Klincke (franziska.klincke@bio.ku.dk)  
   Joseph Nesme (joseph.nesme@bio.ku.dk)  
   Søren Johannes Sørensen (sjs@bio.ku.dk)
