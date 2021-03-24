#######################################################################################################
# !/usr/bin/env python3
# -*- coding:utf-8 -*-
# pipeline for plasmid and mobile
# 2019-03-06
# wanli he
###########################################################################################################

ruleorder: cutting_circular_plasmid>cdhit_nucler_circular_plasmid>circular_index_bam_plasmid>\
           circular_get_bam_file_plasmid>circular_get_sort_file_plasmid>circular_reads_mapping>\
           filter_dectect_circular_contig>circular_relative_contig_abundance>\
           circular_paste_relative_contig_abundance>circular_split_contigs>\
           circular_run_mob_typer>circular_clean_mob_typer_file>circular_makefaa

rule cutting_circular_plasmid:
    input:
        f = expand("circular/plasmid/{sample}_verify_plasmid_circular.fasta", sample=config["samples"])
    output:
        f =temp("circular_non_redundant_contig/all_plasmid_contigs.fa")
    shell:
        "cat {input.f} >>{output.f}"


rule cdhit_nucler_circular_plasmid:
    input:
        f1 = "circular_non_redundant_contig/all_plasmid_contigs.fa"
    output:
        f1 = "circular_non_redundant_contig/circular_non_redundant_contigs.fa"
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        out = "log/circular_non_redundant_contig/cdhit_plasmid.out",
        err = "log/circular_non_redundant_contig/cdhit_plasmid.err"
    params:
        cd = config['cdhit-est_path'],
        c = config['c_cd-hit_c'],
        aS = config['c_cd-hit_aS'],
        M = config['c_cd-hit_M']
    shell:
        "{params.cd} -i {input.f1} " \
                    "-c {params.c} " \
                    "-M {params.M} " \
                    "-aS {params.aS} " \
                    "-o {output.f1} " \
                    "-T {threads}" \
                    " 2>{log.err} >{log.out}"


rule circular_index_bam_plasmid:
    input:
        f = "circular_non_redundant_contig/circular_non_redundant_contigs.fa"
    output:
        f1= temp("circular_non_redundant_contig/circular_non_redundant_contigs.fa.pac"),
        f2= temp("circular_non_redundant_contig/circular_non_redundant_contigs.fa.amb"),
        f3= temp("circular_non_redundant_contig/circular_non_redundant_contigs.fa.bwt"),
        f4= temp("circular_non_redundant_contig/circular_non_redundant_contigs.fa.sa"),
        f5= temp("circular_non_redundant_contig/circular_non_redundant_contigs.fa.ann")
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        err = "log/circular_index_bam_plasmid/index.err",
        out = "log/circular_index_bam_plasmid/index.out"
    shell:
        "bwa index {input.f} ;"
        "rm -rf circular_non_redundant_contig/circular_non_redundant_contigs.fa.clstr"


rule circular_get_bam_file_plasmid:
    input:
        f1= "circular_non_redundant_contig/circular_non_redundant_contigs.fa",
        f2= "qc_reads/{sample}_qc_1.fastq.gz",
        f3= "qc_reads/{sample}_qc_2.fastq.gz",
        f4= "circular_non_redundant_contig/circular_non_redundant_contigs.fa.pac",
        f5= "circular_non_redundant_contig/circular_non_redundant_contigs.fa.amb",
        f6= "circular_non_redundant_contig/circular_non_redundant_contigs.fa.bwt",
        f7= "circular_non_redundant_contig/circular_non_redundant_contigs.fa.sa",
        f8= "circular_non_redundant_contig/circular_non_redundant_contigs.fa.ann"
    output:
        f= temp("circular_non_redundant_contig/{sample}_plasmids.bam")
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        err = "log/circular_get_bam_file_plasmid/{sample}.err",
        out = "log/circular_get_bam_file_plasmid/{sample}.out"
    shell:
        "bwa mem -M {input.f1} {input.f2} {input.f3} -t {threads} | samtools view -Sb - >{output.f}"
        #"bwa mem -M {input.f1} {input.f2} {input.f3} -t {threads}" \
        #"| samtools view -Sb -@ {threads} - >{output.f}"


rule circular_get_sort_file_plasmid:
    input:
        f="circular_non_redundant_contig/{sample}_plasmids.bam"
    output:
        f=temp("circular_non_redundant_contig/{sample}_plasmids_sort.bam")
    threads: config['threads']
    params:
        q = config["bwa_q"]
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    shell:   #python {BWAPAIR}|
        # "samtools view -@ {threads} -h {input.f}|" \
        # "samtools view -Su -q{params.q} -@ {threads} -|samtools sort -@ {threads} -O BAM -o {output.f} -  "
        "samtools view -h {input.f}|" \
        "samtools view -@ {threads} -Su -q{params.q} - |samtools sort  -@ {threads} -O BAM -o {output.f} -"

rule circular_reads_mapping:
    input:
        f="circular_non_redundant_contig/{sample}_plasmids_sort.bam"
    output:
        f=temp("circular_non_redundant_contig/mapping/{sample}_bp_mapping.txt")
    threads: 1
    conda:
        "%s/bedtools.yaml" % CONDAENV
    shell:
        "genomeCoverageBed -ibam {input.f} -d > {output.f}"

rule filter_dectect_circular_contig:
    input:
        f = "circular_non_redundant_contig/mapping/{sample}_bp_mapping.txt",
        f2 = "circular_non_redundant_contig/circular_non_redundant_contigs.fa"
    output:
        f = temp("circular_non_redundant_contig/{sample}_coverage.txt")
    run:
        dict_len={}
        with open(input.f2,"r") as infile2:
            for line in infile2:
                if line.startswith(">"):
                    k = line.strip().split()[0][1:]
                    dict_len[k] = ""
                else:
                    dict_len[k] += line.strip()

        # dict_c = defaultdict(float)
        dict_det = defaultdict(float)
        with open(input.f,"r") as infile:
            for line in infile:
                lst = line.strip().split()
                # dict_c[lst[0]] += float(lst[2])

                if float(lst[2]) != 0:
                    dict_det[lst[0]] += 1
                else:
                    dict_det[lst[0]] += 0

        with open(output.f,"w") as outfile:
            for k in dict_det.keys():
                t = dict_det[k]/float(len(dict_len[k]))
                if t > float(config["contig_detection"]):
                    st = "{}\t{}\n".format(k,t)
                    outfile.write(st)
                else:
                    st = "{}\t{}\n".format(k,0)
                    outfile.write(st)


rule circular_relative_contig_abundance:
    input:
        f = "circular_non_redundant_contig/{sample}_coverage.txt"
    output:
        f = temp("circular_non_redundant_contig/{sample}_relative_abundance.txt")
    threads: 1
    run:
        a = 0
        with open(input.f,"r") as infile:
            for line in infile:
                lst = line.strip().split()
                a += float(lst[1])

        with open(input.f,"r") as infile:
            with open(output.f,"w") as outfile:
                for line in infile:
                    lst = line.strip().split()
                    if a != 0:
                        s= float(lst[1])/float(a)
                        st = lst[0] +"\t" + str(s) +"\n"
                        outfile.write(st)
                    else:
                        st = lst[0] +"\t" +"0\n"
                        outfile.write(st)


rule circular_paste_relative_contig_abundance:
    input:
        f = expand("circular_non_redundant_contig/{sample}_relative_abundance.txt",sample=config["samples"])
    output:
        f = "circular_non_redundant_contig/relative_contig_abundance/all_samples_contig_relative_abundance.txt"
    threads: 1
    run:
        if os.path.exists(output.f):
            pass
        else:
            dict = defaultdict(str)
            header = "contig_ID"

            for file in input.f:
                sample_ID = os.path.basename(file)[:-len("_relative_abundance.txt")]
                head = "\t" + sample_ID
                header += head
                with open(file, "r") as infile:
                    for line in infile:
                        lst = line.strip().split()
                        lt = "\t" + lst[1]
                        dict[lst[0]] += lt

            with open(output.f, "w") as outfile:
                header += "\n"
                outfile.write(header)
                for k in dict.keys():
                    st = "{}{}\n".format(k, dict[k])
                    outfile.write(st)

# classify
rule circular_split_contigs:
    input:
        f1="circular_non_redundant_contig/circular_non_redundant_contigs.fa"
    output:
        f2=temp(directory("circular_plasmid_contig_split"))
    threads: config["threads"]
    run:
        if not os.path.exists(output.f2):
            os.mkdir(output.f2)
        with open(input.f1,mode='r') as infile:
            dir = defaultdict(str)
            for line in infile:
                if line.startswith(">"):
                    k = line
                    dir[k] = ""
                else:
                    dir[k] += line

            for k, v in dir.items():
                res = k[1:].split()[0] + "_split.fa"
                a = os.path.join(output.f2, res)
                with open(a,mode="w") as outfile:
                    outfile.write(k)
                    outfile.write(v)


rule circular_run_mob_typer:
    input:
        f1="circular_plasmid_contig_split"
    output:
        f2=temp(directory("circular_non_redundant_contig/all_result_mob_typer"))
    threads: config["threads_mob-typer"]
    conda:
        "%s/mobtyper.yaml" % CONDAENV
    log:
        err="log/circular_run_classify_plasmid/run_mob_typer.err",
        out="log/circular_run_classify_plasmid/run_mob_typer.out"
    script:
        "../scripts/run_mob_typer.py"



rule circular_clean_mob_typer_file:
    input:
        f1="circular_non_redundant_contig/all_result_mob_typer"
    output:
        f2="circular_non_redundant_contig/plasmid_classify/plasmid_classify.txt"
    threads: config["threads"]
    run:
        path=input.f1
        files=glob.glob(path+"/*_mob_typer_result")
        lst=[]
        header=""
        for file in files:
            son_path = file+"/*"
            for son_file in glob.glob(son_path):
                with open(son_file,mode="r") as infile:
                    header = infile.readline()
                    st=infile.readline()
                    lst.append(st)

        with open(output.f2,mode="w") as outfile:
            outfile.write(header)
            for ele in lst:
                outfile.write(ele)


#co-exist
rule circular_makefaa:
    input:
        f="circular_non_redundant_contig/circular_non_redundant_contigs.fa"
    output:
        f1 = temp("circular_non_redundant_contig/all_genecalling_protein.faa"),
        f2 = temp("circular_non_redundant_contig/all_genecalling_nucl.fa")
    threads: config["threads"]
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        out = "log/linear_co-exist_circular/all_genecalling.out",
        err = "log/linear_co-exist_circular/all_genecalling.err"
    params:
        p = config["c_prodigal_p"]
    shell:
        "prodigal -p {params.p} -i {input.f} -a {output.f1} -d {output.f2} 2>{log.err} >{log.out}"


rule circular_plasmid_with_MGEs:
    input:
        f = "circular_non_redundant_contig/all_genecalling_protein.faa"
    output:
        f = "circular_non_redundant_contig/co-exist/all_plasmids_with_MEGs.txt"
    threads: config["threads"]
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    params:
        mgesdb = config["MGEs_database"]
    log:
        out = "log/co-exist/linear_co-exist-mges.out",
        err = "log/co-exist/linear_co-exist-mges.err"
    shell:
        "hmmsearch --tblout {output.f} {params.mgesdb} {input.f} 2>{log.err} >{log.out}"


rule circular_contig_annotation_ARGs:
    input:
        f = "circular_non_redundant_contig/circular_non_redundant_contigs.fa"
    output:
        f = directory("circular_non_redundant_contig/co-exist/annotation_ARGs")
    threads: config["threads"]
    conda:
        "%s/rgi.yaml" % CONDAENV
    log:
        out = "log/co-exist/circular_ARGs.out",
        err = "log/co-exist/circular_ARGs.err"
    params:
        d = config['c_rgi_d'],
        it = config['c_rgi_input_type'],
        a = config['c_rgi_a'],
        mode = config["c_rgi_mode"]
    shell:
        "mkdir {output.f};" \
        "rgi main -i {input.f} " \
                "--output_file {output.f} " \
                "-d {params.d} " \
                "--input_type {params.it} " \
                "-n {threads} " \
                "-a {params.a} " \
                "{params.mode} --clean " \
                "2>{log.err} >{log.out};" \
        "mv {output.f}.json {output.f}.txt {output}"


rule circular_contig_annotation_vf:
    input:
        f = "circular_non_redundant_contig/all_genecalling_protein.faa"
    output:
        f = "circular_non_redundant_contig/co-exist/annotation_vf/annotation_vf.txt"
    threads: config["threads"]
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    params:
        db = config["vfdb_database"],
        outfmt = config["c_diamond_outfmt_vf"],
        evalue = config["c_diamond_evalue_vf"],
        id = config["c_diamond_id_vf"],
        cover = config["c_diamond_cover_vf"],
        md = config["c_diamond_blast_module_vf"]
    log:
        out = "log/co-exist/circular_vf.out",
        err = "log/co-exist/circular_vf.err"
    shell:
        "diamond {params.md} " \
                "--db {params.db} " \
                "--query {input.f} " \
                "--out {output.f} " \
                "--outfmt {params.outfmt} " \
                "--evalue {params.evalue} " \
                "--id {params.id} " \
                "--query-cover {params.cover} " \
                "--threads {threads} " \
                "2>{log.err} >{log.out}"

rule circular_contig_annotation_BacMet2_gene_database:
    input:
        f = "circular_non_redundant_contig/all_genecalling_protein.faa"
    output:
        f = "circular_non_redundant_contig/co-exist/annotation_BacMet2/annotation_BacMet2.txt"
    threads: config["threads"]
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    params:
        db = config["BacMet2_database"],
        outfmt = config["c_diamond_outfmt_BacMet2"],
        evalue = config["c_diamond_evalue_BacMet2"],
        id = config["c_diamond_id_BacMet2"],
        cover = config["c_diamond_cover_BacMet2"],
        md = config["c_diamond_blast_module_BacMet2"]
    log:
        out = "log/co-exist/circular_BacMet2.out",
        err = "log/co-exist/circular_BacMet2.err"
    shell:
        "diamond {params.md} " \
                "--db {params.db} " \
                "--query {input.f} " \
                "--out {output.f} " \
                "--outfmt {params.outfmt} " \
                "--evalue {params.evalue} " \
                "--id {params.id} " \
                "--query-cover {params.cover} " \
                "--threads {threads} " \
                "2>{log.err} >{log.out}"

#gene synteny
rule circular_gene_synteny:
    input:
        f = "circular_non_redundant_contig/all_genecalling_nucl.fa"
    output:
        f = "circular_non_redundant_contig/gene_synteny.txt"
    run:
        dict = defaultdict(str)
        with open(input.f,"r") as infile:
            for line in infile:
                if line.startswith(">"):
                    # if lst[8].split(";")[1] == "partial=00":
                    lst = line.strip().split()
                    cid = lst[0][1:-2]
                    gid = lst[0].rsplit("_",1)[1]
                    start = lst[2]
                    end = lst[4]
                    dir = lst[6]
                    v = "({},{},{},{})".format(gid,start,end,dir)
                    dict[cid] += v+"-"

        with open(output.f,"w") as outfile:
            for k in dict.keys():
                st = "{}\t{}\n".format(k,dict[k].strip("_"))
                outfile.write(st)