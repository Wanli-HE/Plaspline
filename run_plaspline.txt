#--------------------------- create a conda environment and insatll dependencies:
     conda create -n plaspline -c conda-forge -c bioconda python=3.9 biopython=1.78 snakemake=5.25.0 click=7.1.2 ruamel.yaml=0.16.12

#---------------------------Clone package
    git clone https://github.com/Wanli-HE/Plaspline.git 

#-------------------------------------------downloading tools and database
     python <~/Plaspline/bin/run_main.py> downloading  -j 5
 
     "-j 5" means 5 job parallel.
     <~/Plaspline/bin/run_main.py> is must the absoulate path of file "run_main.py"

#---------------------------------------------------prepare config.yaml and sample.json file.
     python  <~/Plaspline/bin/run_main.py> preprocessing --path <your reads folder> --adapter <your adapter file> --phix  <your phix file> --assembler <defult:'spades'> --skip_qc <defult:"False">

Notice: 
1. In <your reads folder>: plaspline can only recognize the raw reads (forward and resvers) file which suffix is "_R1","_r1"  or "1" (only show forward in example).  
      *****must be gzip file. otherwise in "quality control step" it will raise error.

2. If you want to do the quality control step,  <your adapter file> and <your phix file> must be input;    if your input is qc reads, not necessary do the  quality control step,  you don`t need input <your adapter file> and <your phix file>, but the  "--skip_qc" must be "True".

3. the assembler defult is spades, you can also choose megahit as assembler, by input the params "--assembler megahit".

for more information, "python Plaspline/bin/run_main.py preprocessing -h" 


# ------------------------------------------ working
1. before true runing plaspline, it is better to check whether there is any command line errors, by dry run (" -n ").
       python <~/Plaspline/bin/run_main.py> working  "step" -j 5  -n 

       the "-j 5" means 5 job parallel.   here notice the total threads that you are using, in the "config.yaml" file records the max threads of each rule ("threads: 8" in config.yaml). so the total threads = (threads)*(job) = 5*8 =40. so make sure that you have enough threads in your system.

2. True runing:
        python <~/Plaspline/bin/run_main.py> working  "step" -j 5 

3. if you are working on cluster system. there are two to run the command:
        a:   python <~/Plaspline/bin/run_main.py> working "step"  -j 5  --cluster 'qsub -t 40 -l nodes=1'  
        b:   python <~/Plaspline/bin/run_main.py> working "step"  -j 5  --profile cluster --latency-wait 60
        
        the method b is utilized Altas (https://metagenome-atlas.readthedocs.io/en/latest/usage/getting_started.html#usage), so user can excaute Plaspline in cluster system according to "Execue Atlas / Cluster execution" module.

Runing plaspline in each steps:
1.quality control step
        python <~/Plaspline/bin/run_main.py> working "qc"  -j 5
	
   Notice:
        the raw reads file need be "gzip" file, otherwise is will raise error in here and give a hint "check you raw and qc reads fastq file, whether is gzip file!"
        and also it is very import that confirm and modify the params in config.yaml file. 
	#------ rule adapter_moving params
	     adapter_bbduk_r: r
	     adapter_bbduk_k: 23
	     adapter_bbduk_mink: 11
	     adapter_bbduk_mem: 5g
	     adapter_file: ""
	#------ rule phix_moving params
	     phix_k: 31
	     phix_hdist: 1
	     phix_mem: 5g
	     phix_file: ""
	#------ rule quality_window_trimming params
	     filter_outgzip: -g
	     filter_t: sanger
	     filter_l: 50
	     filter_m: pe

2. assembly step
       python <~/Plaspline/bin/run_main.py> working "assembly"  -j 5

       Notice:
        it is import that the max k-mer when doing assembly. because in the next step "circular" part which need the max k-mer information.
        if the assebler is "spades", the params are in config.yaml file. 
	spades_filter_contig: 500
	spades_k: auto                   #just change "auto" to a series number for assembly. 21,33,55 for example
	spades_m: 5000
	spades_preset: meta

              **** if the spades_k is auto you need to check the max k-mer in the assembly folder of each sample in which records the max k-mer

        if the assembler is "megahit". you need to confirm the params in config.yaml file. 
	megahit_core_len: 121
	megahit_k_max: 121
	megahit_k_min: 21
	megahit_k_step: 20
	megahit_low_local_ratio: 0.2
	megahit_m: 0.95
	megahit_merge_level: 20,0.98
	megahit_min_count: 2
	megahit_preset: default
	megahit_prune_level: 2

3. circular step
       before runing this step, confim and modify config.yaml of tools params:
                   #metaplasmidspades
	     c_spades_k: auto
	     c_spades_m: 5000
	     c_sapdes_preset1: meta
	     c_sapdes_preset2: plasmid

                   #scapp:
	     scapp_k: 121           #******* here is the max assembly k-mer, must be confirm and modify  
	     scapp_gh: False
        
        then runing the command:
                   python <~/Plaspline/bin/run_main.py> working "circular"  -j 5

4.isolation linear plasmid contigs from all assembled contigs
       python <~/Plaspline/bin/run_main.py> working "isolation"  -j 5
       
       Notice:     confim and modify config.yaml of tools params
	#plasflow
	plasflow_threshold: 0.7

	#plton
	platon_mode: accuracy
        	  	 # sensitivity   RDS>0.95,
        		   # accuracy      RDS & characterization heuristics, highest accuracy -----default
          		 # specificity   RDS>0.999

5. non-redundant circular plasmid contig set step
      python <~/Plaspline/bin/run_main.py> working "contig_circular"  -j 5

      Notice:     confim and modify config.yaml of tools params
            # cluster all plasmid contigs	
	circular_contig_mmseqs2_id: 0.9
	circular_contig_mmseqs2_c: 0.95
	circular_contig_mmseqs2_mode: 0

           #functional_annotation
	c_functional_annotation_m: diamond

           #mges gene annotation
	MGEs_database: /mibi/Wanli/MEGs.database/hgtidb_pfams.hmm

           #ARGs gene annotation
	c_rgi_d: plasmid   #[wgs, plasmid, chromosome, NA]
	c_rgi_input_type: contig   #[contig, protein]
	c_rgi_a: DIAMOND
	c_rgi_mode:  "" #[None, --include_loose, --exclude_nudge, --include_loose --exclude_nudge ]
               		 #https://github.com/arpcard/rgi   ---- for more details

            #c_vfdb_gene_database: #path
	c_diamond_outfmt_vf: 6
	c_diamond_evalue_vf: 1e-5
	c_diamond_id_vf: 70
	c_diamond_cover_vf: 0.7
	c_diamond_blast_module_vf: blastp

             #BacMet2 gene annotation
	c_diamond_outfmt_BacMet2: 6
	c_diamond_evalue_BacMet2: 1e-5
	c_diamond_id_BacMet2: 70
	c_diamond_cover_BacMet2: 0.7
	c_diamond_blast_module_BacMet2: blastp


6. contig circular classify step
      python <~/Plaspline/bin/run_main.py> working "contig_circular_classify"  -j 5
     
      Notice:     confim and modify config.yaml of tools params
               #classify
	threads_mob-typer: 24

7. non-redundantlinear plasmid contig set step
      python <~/Plaspline/bin/run_main.py> working "contig_linear"  -j 5
	
      Notice:     confim and modify config.yaml of tools params
               # cluster
	linear_contig_mmseqs2_id: 0.9
	linear_contig_mmseqs2_c: 0.95
	linear_contig_mmseqs2_mode: 1
               
            other params in here are same with  " non-redundant circular plasmid contig set step " (below)
           #functional_annotation
	c_functional_annotation_m: diamond

           #mges gene annotation
	MGEs_database: /mibi/Wanli/MEGs.database/hgtidb_pfams.hmm

           #ARGs gene annotation
	c_rgi_d: plasmid   #[wgs, plasmid, chromosome, NA]
	c_rgi_input_type: contig   #[contig, protein]
	c_rgi_a: DIAMOND
	c_rgi_mode:  "" #[None, --include_loose, --exclude_nudge, --include_loose --exclude_nudge ]
               		 #https://github.com/arpcard/rgi   ---- for more details

            #c_vfdb_gene_database: #path
	c_diamond_outfmt_vf: 6
	c_diamond_evalue_vf: 1e-5
	c_diamond_id_vf: 70
	c_diamond_cover_vf: 0.7
	c_diamond_blast_module_vf: blastp

             #BacMet2 gene annotation
	c_diamond_outfmt_BacMet2: 6
	c_diamond_evalue_BacMet2: 1e-5
	c_diamond_id_BacMet2: 70
	c_diamond_cover_BacMet2: 0.7
	c_diamond_blast_module_BacMet2: blastp
         
8. circular plasmid gene analysis step
      python <~/Plaspline/bin/run_main.py> working "gene_circular"  -j 5

     Notice:     confim and modify config.yaml of tools params
            #reads filter
	circular_msamtools_gene_p: 90      #the percent identity of alignment is >=90%
	circular_msamtools_gene_l: 80      #the read is at least 80bp long
	circular_msamtools_gene_z: 80      #at least 80% of the read is aligned

            #gene abundance
	g_prodigal_p: meta

	g_cd-hit_c: 0.9
	g_cd-hit_aS: 0.95
	g_cd-hit_M: 20000

	bwa_q: 15

            #gene abundance or relative abundance
	circular_gene_abundance_format: ab     #rel is relative abundance; ab is abundance

           #functional_annotation
	g_functional_annotation_m: diamond	
                  functional_threads: 24
           #ARGs gene annotation
	g_rgi_d: plasmid   #[wgs, plasmid, chromosome, NA]
	g_rgi_input_type: contig   #[contig, protein]
	g_rgi_a: DIAMOND
	g_rgi_mode: "" #[None, --include_loose, --exclude_nudge, --include_loose --exclude_nudge ]
                		#https://github.com/arpcard/rgi   ---- for more details

          #c_vfdb_gene_database: #path
	g_diamond_outfmt_vf: 6
	g_diamond_evalue_vf: 1e-5
	g_diamond_id_vf: 70
	g_diamond_cover_vf: 0.7
	g_diamond_blast_module_vf: blastp


         #BacMet2 gene annotation
	g_diamond_outfmt_BacMet2: 6
	g_diamond_evalue_BacMet2: 1e-5
	g_diamond_id_BacMet2: 70
	g_diamond_cover_BacMet2: 0.7
                   g_diamond_blast_module_BacMet2: blastp

	
9.  linear plasmid gene analysis step
      python <~/Plaspline/bin/run_main.py> working "gene_linear"  -j 5

      Notice:     confim and modify config.yaml of tools params

            #reads filter
	linear_msamtools_gene_p: 90      #the percent identity of alignment is >=90%
	linear_msamtools_gene_l: 80      #the read is at least 80bp long
	linear_msamtools_gene_z: 80      #at least 80% of the read is aligned

            #gene abundance
	g_prodigal_p: meta

	g_cd-hit_c: 0.9
	g_cd-hit_aS: 0.95
	g_cd-hit_M: 20000

	bwa_q: 15

            #gene abundance or relative abundance
	circular_gene_abundance_format: ab     #rel is relative abundance; ab is abundance

           #functional_annotation
	g_functional_annotation_m: diamond	
                  functional_threads: 24
           #ARGs gene annotation
	g_rgi_d: plasmid   #[wgs, plasmid, chromosome, NA]
	g_rgi_input_type: contig   #[contig, protein]
	g_rgi_a: DIAMOND
	g_rgi_mode: "" #[None, --include_loose, --exclude_nudge, --include_loose --exclude_nudge ]
                		#https://github.com/arpcard/rgi   ---- for more details

          #c_vfdb_gene_database: #path
	g_diamond_outfmt_vf: 6
	g_diamond_evalue_vf: 1e-5
	g_diamond_id_vf: 70
	g_diamond_cover_vf: 0.7
	g_diamond_blast_module_vf: blastp


         #BacMet2 gene annotation
	g_diamond_outfmt_BacMet2: 6
	g_diamond_evalue_BacMet2: 1e-5
	g_diamond_id_BacMet2: 70
	g_diamond_cover_BacMet2: 0.7
                   g_diamond_blast_module_BacMet2: blastp




 	