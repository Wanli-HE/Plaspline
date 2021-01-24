# plasmid-pipeline

# Introduction
Nowadays, Numerous tools are developed based on different strategies for detecting plasmid from metagenome. Apparently, all of them have merits, but the demerits are also non-negligible. Thus, how to estimate them and how to combine the advantages of these tools are necessary. Besides, it is lack of a proven, systematic, and comprehensive workflow for theses multiple preprocessing and analytical steps, which are not only including detecting plasmid, but also including significant downstream analysis of plasmid.   
  

we are benchmarking these tools and strategies, building up *Plaspline* based on the benchmarked tools, which is trying to combine the advantages of tools for a better result of detecting plasmid. Moreover, multiple downstream analysis is also added into this workflow includes quality control, assembly, circular/linear plasmids isolation, circular/linear plasmids genome verification and classification, and plasmids-relative-genes analysis, which aims to get a comprehensive analysis of plasmid, both in gene and plasmid community level. 


# Requirements
ruamel.yaml  
click  
snakemake  
biopython

# Install
Clone package:

    > git clone https://github.com/Wanli-HE/plasmid-pipeline.git  
    > cd plasmid-pipeline

all scripts are under the folder

# Usage
1. all databases and dependent tools are downloaded and installed automatically with default version by script:  
    
         > python ~/plasmid-pipeline/bin/run_main.py downloading 
         
     
   before running this command, you need to install anaconda or miniconda. If you haven’t done it already you need to configure conda with the bioconda-channel and the conda-  forge channel. This are sources for packages beyond the default one:

        > conda config --add channels defaults  
        > conda config --add channels bioconda  
        > conda config --add channels conda-forge

   the version of databastes and tools is recorded in "dependdent-tools-database.txt" file,  and these tools are installed by two ways, one is using conda yaml file, another is using "wget" command line (all is recored in "dependdent-tools-database.txt").

   NOTICE: if not necessary, suggesting you do not change the version of the tools, because we do not know what had changed in the new version of the tool and whether it will have any conflict with Plaspline.  
     
   the tools and databases which installed by "wget+web-downloading-address", you can using "python ~/plasmid-pipeline/bin/run_main.py downloading -h" to check, and also can change the peremters to choose the version that you need.  
     
   the tools and databases which installed by conda YAML file, in ~/plasmid-pipeline/envs, you can also change the version from there.  
   
2. configurating your config.yaml and sample.josn files, using command:  
   
       > python ~/plasmid-pipeline/bin/run_main.py preprocessing    

   the config.yaml and sample.josn files will generate automatically in current folder. *the later steps and commands will excaute only in this folder, otherwise it will raise error which "can not find config.yaml and sample.json".   
     
   once config.yaml are generated, user can configurate it manually, if any params which you can not understand, you can chenck "the default_config.yaml" file which is in the ~/plasmid-pipeline/conf folder. before you running this script, you need to configurate your own config.yaml.   
     
   sample.json file is created also in this command, but onething is important is your reads file suffix, plaspline can only recognize "_R1","_r1","_1". because i only meet these three suffixs -_-, so before you run this command, making sure your reads suffix can be recognize by Plaspline. if not, you can changing it by yourself or contact with us, we will upgrade it.  
   
3. running the pipeline
    using the script:  
    
       > python ~/plasmid-pipeline/bin/run_main.py working 
       
    and before running, REMEMBER!!!!! checking your congig.yaml file!!!!!!


# Contact

## author:   
   Wanli HE (wanli.he@bio.ku.dk)

## cooperator:  
   Franziska Klincke (franziska.klincke@bio.ku.dk)  
   Joseph Nesme (joseph.nesme@bio.ku.dk)  
   Søren Johannes Sørensen (sjs@bio.ku.dk)
