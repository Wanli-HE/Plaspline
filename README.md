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


# Usage
1. all databases and dependent tools are downloaded and installed automatically with default version by script:  
    
         > python ~/Plaspline/bin/run_main.py downloading   
         
     
   before running this command, you need to install anaconda or miniconda. If you haven’t done it already you need to configure conda with the bioconda-channel and the conda-  forge channel. This are sources for packages beyond the default one:

        > conda config --add channels defaults  
        > conda config --add channels bioconda  
        > conda config --add channels conda-forge

   these tools are installed by two ways, one is using conda yaml file, another is using "wget + download-address" command line.

   NOTICE: if not necessary, suggesting you do not change the version of the tools, because we do not know what had changed in the new version of the tool and whether it will have any conflict with Plaspline.  
     
   the tools and databases which installed by "wget+web-downloading-address", you can using "python ~/plasmid-pipeline/bin/run_main.py downloading -h" to check, and also can change the peremters to choose the version that you need.  
     
   the tools and databases which installed by conda YAML file, in ~/plasmid-pipeline/envs, you can also change the version from there.  
   
2. configurating your config.yaml and sample.josn files, using command:  
   
       > python ~/Plaspline/bin/run_main.py preprocessing    

   the config.yaml and sample.josn files will generate automatically in current folder. *the later steps and commands will excaute only in this folder, otherwise it will raise error which "can not find config.yaml and sample.json".   
     
   once config.yaml are generated, user can configurate it manually, if any params which you can not understand, you can chenck "the default_config.yaml" file which is in the ~/plasmid-pipeline/conf folder. before you running this script, you need to configurate your own config.yaml.  
   
   NOTICE: if the assembler is SPAdes, one parameter you must configurate is the Max assembly k-mer, because it is relative with plasmid detection when SCAPP are running in the later step. the default "-k" of SPAdes is "auto", but one parameter of SCAPP is "-k" which is the Max assembly k-mer. Thus, when you excuate Plaspline, one way is configurate the assembler k-mer (for example changing "spades_k: auto" to "spades_k: [21,33,55]"), and then configurate "scapp_k: 55", 55 is the max assembly k-mer. if you do not know which assembler k-mer you should choose, another way is run the Plaspline dividually, you do not need configurate "spades_K", it is still "auto", but first run assembly and check the assembler max k-mer, and then configurate "scapp_k" to max k-mer.    
     
   sample.json file is created also in this command, but onething is important is your reads file suffix, plaspline can only recognize "_R1","_r1","_1". because i only meet these three suffixs -_-, so before you run this command, making sure your reads suffix can be recognize by Plaspline. if not, you can changing it by yourself or contact with us, we will upgrade it.  
   
3. running the pipeline
    using the command:  
    
       > python ~/Plaspline/bin/run_main.py working 
       
    and before running, REMEMBER!!!!! checking your congig.yaml file!!!!!!  
    
    running all step:  
    
       > python  ~/Plaspline/bin/run_main.py working "all" -j 25
       
    running dividually:  
    
       > python  ~/Plaspline/bin/run_main.py working "assembly" -j 25
       
       (more parameter you can check by python  ~/Plaspline/bin/run_main.py working -h)  
    
    running in cluster system:  
    
       > python  ~/Plaspline/bin/run_main.py working "circular" -j 25  --profile cluster --latency-wait 60  
       
      for running in cluster system, we are utilized  Altas (https://metagenome-atlas.readthedocs.io/en/latest/usage/getting_started.html#usage), so user can excaute Plaspline in cluster system according to "Execue Atlas / Cluster execution" module.  
      
      Notice: "--latency-wait 60" is necessary, otherwise it will raise error "MissingInputException..." (this is a bug of snakemake)
    

# Contact

## author:   
   Wanli HE (wanli.he@bio.ku.dk)

## cooperator:  
   Franziska Klincke (franziska.klincke@bio.ku.dk)  
   Joseph Nesme (joseph.nesme@bio.ku.dk)  
   Søren Johannes Sørensen (sjs@bio.ku.dk)
