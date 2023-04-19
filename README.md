# Plaspline, a pipeline for plasmid analysis form metagenome

# Introduction
*Plaspline* is a snakemake based workflow, which aims to get a comprehensive analysis of plasmid, both in gene and plasmid community level from short-gun sequences. The main steps are quality control, assembly, circular/linear plasmids isolation, circular/linear plasmids genome verification and classification, and plasmids-relative-genes analysis.


# Requirements
ruamel.yaml=0.16.12  
click=7.1.2  
snakemake=5.25.0  
biopython=1.78  
python=3.9

# Install  
create a conda environment and insatll dependencies:  
```
> conda create -n plaspline -c conda-forge -c bioconda python=3.9 biopython=1.78 snakemake=7.22.0 click=7.1.2 ruamel.yaml=0.16.12
```

Clone package:  
```
> git clone https://github.com/Wanli-HE/Plaspline.git 
```  
all scripts are under the 'Plaspline' folder


# Downlanding dependent databases and tools:

*** make sure that you have enough space in your folder for it (80G maybe!!!!)

downloading the databases and tools
```
> python Plaspline/bin/run_main.py downloading       
``` 
*** Notice: when downloading plasforest database, it requires input an e-mail address in your terminal, to get full database or input "no"!!!


# Usage
## Prepare config.yaml and sample.json file.

```
> python  Plaspline/bin/run_main.py preprocessing --path Plaspline/example/reads --adapter Plaspline/example/adapters.fa.gz --phix Plaspline/example/phix174_ill.ref.fa.gz
```
Notice: 
1. In your reads folder: plaspline can only recognize the raw reads (forward and resvers) file which suffix is "_R1" or "_r1" or "_1" (only show forward in example).  
      ***must be gzip file. otherwise in "quality control step" it will raise error.

2. If you want not to do the quality control step, but still have to give the "--adapter" and "--phix" input, (a fake path is ok,  for example ./ ). don`t forget change "--skip_qc " to "True", which means jump quality control step.  

3. the assembler defult is spades, you can also choose megahit as assembler, by input the params "--assembler megahit".
     
4. in command line "Plaspline/bin/run_main.py" is relative path, if you are not running the command with same folder of Plaspline, you need to input absolute path (abs.path/Plaspline/bin/run_main.py)

for more information, "python ~Plaspline/bin/run_main.py preprocessing -h" 
  
## Working
### Dry running 
Before true runing plaspline, it is better to check whether there is any command line errors, by dry run (" -n ").

```
  > python Plaspline/bin/run_main.py working  <"step"> -j 3  -n 
```
Notice:
  1. before running, make sure the "cofig.yaml" and "samples.json" are under the floder which you are running the command.
  2. **** choosing jobs (-j):
     - in config.yaml file, it records {threads} which is each rules threads, so in "-j (pnum)", this pnum <= sample_number * {threads}      
     - running in cluster, -j {pnum};  this pnum <= sample number
          
### True running:
```
> python Plaspline/bin/run_main.py working  "step" -j 3 
```
     
### Runing plaspline all:

*** before running it, you have to make sure that scapp_k == max ker (that you are using in assembly step) in config.yaml file

     
```
> python Plaspline/bin/run_main.py working  "all" -j 3
```
     



     
### Runing plaspline in each steps: (suggestion!!!)

1.quality control step
  ```
  > python Plaspline/bin/run_main.py working "qc"  -j 3
  ```
2. assembly step
  ```
  > python Plaspline/bin/run_main.py working "assembly"  -j 3
  ```
3. circular step
     
  before running this step, check the max k-mer you are using in assembly step, and modifying scapp_k to "max k-mer (assembly)" in config.yaml file.
  ```
  > python Plaspline/bin/run_main.py working "circilar"  -j 3
  ```
4. isolation linear plasmid contigs from all assembled contigs
  ```
  > python Plaspline/bin/run_main.py working "isolation"  -j 3
  ```
5. non-redundant circular plasmid contig set step
  ```
  > python Plaspline/bin/run_main.py working "contig_circular"  -j 3
  ```
6. contig circular classify step
  ```
  > python Plaspline/bin/run_main.py working "contig_circular_classify"  -j 3
  ```
7. non-redundantlinear plasmid contig set step
  ```
  > python Plaspline/bin/run_main.py working "contig_linear"  -j 3
  ```
8. circular plasmid gene analysis step
  ```
  > python Plaspline/bin/run_main.py working "gene_circular"  -j 3
  ```
9. linear plasmid gene analysis step
  ```
  > python Plaspline/bin/run_main.py working "gene_linear"  -j 3
  ```


     
     
# Contact

## author:   
   Wanli HE (wanli.he@bio.ku.dk)

## cooperator:  
   Franziska Klincke (franziska.klincke@bio.ku.dk)  
   Joseph Nesme (joseph.nesme@bio.ku.dk)  
   Søren Johannes Sørensen (sjs@bio.ku.dk)
