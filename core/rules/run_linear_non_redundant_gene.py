#######################################################################################################
# !/usr/bin/env python3
# -*- coding:utf-8 -*-
# pipeline for plasmid and mobile
# 2021-1-21
# wanli he
###########################################################################################################


rule linear_genecalling_plasmid_gene:
    input:
        f = "linear_plasmid_genome/{sample}_predict_plasmid.fa"
    output:
        f = temp("linear_non_redundant_gene/{sample}_nucl_gene.fa")
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    params:
        p = config['g_prodigal_p']
    log:
        out = "log/linear_non_redundant_gene/{sample}_genecalling_genecalling.out",
        err = "log/linear_non_redundant_gene/{sample}_genecalling_genecalling.err"
    script:
        "../scripts/gene_calling.py"


rule linear_cutting_all_plasmid_gene:
    input:
        f = expand("linear_non_redundant_gene/{sample}_nucl_gene.fa", sample=config["samples"])
    output:
        f =temp("linear_non_redundant_gene/all_plasmid_genes.fa")
    shell:
        "cat {input.f} >>{output.f}"


rule rename_genecalling_file_ff:
    input:
        f = "linear_non_redundant_gene/all_plasmid_genes.fa",
    output:
        f = temp("linear_non_redundant_gene/all_contigs_nucl_gene_rename.fa")
    run:
        with open(input.f,"r") as infile:
            with open(output.f,"w") as outfile:
                for line in infile:
                    if line.startswith(">"):
                        lst = line.strip().split("\t",1)
                        la = lst[0].rsplit("_",1)[0]
                        st = la+"\n"
                        outfile.write(st)
                    else:
                        outfile.write(line)


rule non_redundant_nucler_linear_gene:
    input:
        f = "linear_non_redundant_gene/all_contigs_nucl_gene_rename.fa"
    output:
        f = "linear_non_redundant_gene/linear_non_redundant_gene"
    threads: config['threads']
    conda:
        "%s/mmseqs2.yaml" % CONDAENV
    log:
        out = "log/linear_non_redundant_gene/non_redundant_plasmid.out",
        err = "log/linear_non_redundant_gene/non_redundant_plasmid.err"
    params:
        id = config["linear_gene_mmseqs2_id"],
        c =  config["linear_gene_mmseqs2_c"],
        mode =  config["linear_gene_mmseqs2_mode"]
    shell:
        "mmseqs easy-cluster {input.f} {output.f} tmp --min-seq-id {params.id} -c {params.c} " \
        "--cov-mode {params.mode} 2>{log.err} >{log.out};" \
        "touch {output.f}"


rule linear_gene_ab:
    input:
        f1 = "linear_non_redundant_contig/contig_abundance/per_sample/{sample}_existing_contigs_ab_per_sample.txt",
        f3 = "linear_non_redundant_gene/linear_non_redundant_gene",
        f2 = "linear_non_redundant_contig/linear_non_redundant_contigs"
    output:
        f = temp("linear_non_redundant_gene/{sample}_ab_final.txt")
    run:
        dict_rb = {}  # k rep ; v ab   contig
        with open(input.f1, 'r') as infile:
            for line in infile:
                lst = line.strip().split()
                dict_rb[lst[0]] = float(lst[1])

        dict_contig_catalog = {}   # k query ; v rep
        with open(input.f2+"_cluster.tsv",'r') as infile:
            for line in infile:
                lst = line.strip().split()
                dict_contig_catalog[lst[1]] = dict_rb[lst[0]]

        dict_gene_catalog = defaultdict(list)     # k : rep_gene ; v [contig_query1,contig_query2]
        with open(input.f3+'_cluster.tsv','r') as infile:
            for line in infile:
                lst = line.strip().split()
                dict_gene_catalog[lst[0]].append(lst[1].rsplit("_",1)[0])

        with open(output.f,'w') as outfile:
            for k in dict_gene_catalog.keys():
                b = 0.0
                for i in dict_gene_catalog[k]:        # dict_gene_catalog[k] is contig id [contig_query1,contig_query]
                    c = dict_contig_catalog[i]
                    b += float(c)
                st = "{}\t{}\n".format(k,b)
                outfile.write(st)


rule cutting_all_linear_coverage_file:
    input:
        f = expand("linear_non_redundant_gene/{sample}_ab_final.txt",sample=config["samples"])
    output:
        f = "linear_non_redundant_gene/abundance/all_linear_gene_abundance_based_contig_abundance.txt"
    run:
        if os.path.exists(output.f):
            pass
        else:
            dict = defaultdict(str)
            header = "gene_ID"

            for file in input.f:
                sample_ID = os.path.basename(file)[:-len("_ab_final.txt")]
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
                    st = "{}\t{}\n".format(k, dict[k])
                    outfile.write(st)


#gene functional annotation
rule linear_gene_genecalling_gene:     #non-redundant-gene-set
    input:
        f = "linear_non_redundant_gene/linear_non_redundant_gene"
    output:
        f = temp("linear_non_redundant_gene/linear_gene_prodigal_protein_seq__.faa")
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    params:
        p = config["g_prodigal_p"]
    log:
        out ="log/linear_gene_annotation/plasmid_gene_genecalling.out",
        err ="log/linear_gene_annotation/plasmid_gene_genecalling.err"
    shell:
        "prodigal -p {params.p} " \
                "-i {input.f}_rep_seq.fasta " \
                "-a {output.f} " \
                "2>{log.err} >{log.out}"

rule rename_genecalling_file_cc_bb:
    input:
        f = "linear_non_redundant_gene/linear_gene_prodigal_protein_seq__.faa",
    output:
        f = temp("linear_non_redundant_gene/linear_gene_prodigal_protein_seq.faa")
    run:
        with open(input.f,"r") as infile:
            with open(output.f,"w") as outfile:
                for line in infile:
                    if line.startswith(">"):
                        lst = line.strip().split("\t",1)
                        la = lst[0].rsplit("_",1)[0]
                        st = la+"\n"
                        outfile.write(st)
                    else:
                        outfile.write(line)


rule linear_functional_annotation_genes:
    input:
        f = "linear_non_redundant_gene/linear_non_redundant_gene"
    output:
        f = directory("linear_non_redundant_gene/functional_annotation")
    threads: config['functional_threads']
    conda:
        "%s/emapper.yaml" % CONDAENV
    log:
        out = "log/linear_gene_annotation/functional_annotation.out",
        err = "log/linear_gene_annotation/functional_annotation.err"
    params:
        m = config['g_functional_annotation_m'],
        d = E_DB
    script:
        "../scripts/emapper.py"


rule linear_plasmid_gene_with_MGEs:
    input:
        f = "linear_non_redundant_gene/linear_gene_prodigal_protein_seq.faa"
    output:
        f = "linear_non_redundant_gene/MGES-annotation/all_plasmids_with_MEGs.txt"
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

#other gene annotation
rule linear_annotation_ARGs_gene:
    input:
        f = "linear_non_redundant_gene/linear_non_redundant_gene"
    output:
        f = directory("linear_non_redundant_gene/annotation_ARGs")
    threads: config["threads"]
    conda:
        "%s/rgi.yaml" % CONDAENV
    log:
        out = "log/linear_gene_annotation/annotation_ARGs.out",
        err = "log/linear_gene_annotation/annotation_ARGs.err"
    params:
        d = config['g_rgi_d'],
        it = config['g_rgi_input_type'],
        a = config['g_rgi_a'],
        mode = config["g_rgi_mode"]
    shell:
        "mkdir -p {output.f};" \
        "rgi main -i {input.f}_rep_seq.fasta " \
                "--output_file {output.f} " \
                "-d {params.d} " \
                "--input_type {params.it} " \
                "-n {threads} " \
                "-a {params.a} " \
                "{params.mode} --clean  " \
                "2>{log.err} >{log.out};" \
        "mv {output.f}.json {output.f}.txt {output}"


rule linear_annotation_vf_gene:
    input:
        f = "linear_non_redundant_gene/linear_gene_prodigal_protein_seq.faa"
    output:
        f = "linear_non_redundant_gene/annotation_vf/annotation_vf.txt"
    threads: config["threads"]
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    params:
        db = config["vfdb_database"],
        outfmt = config["g_diamond_outfmt_vf"],
        evalue = config["g_diamond_evalue_vf"],
        id = config["g_diamond_id_vf"],
        cover = config["g_diamond_cover_vf"],
        md = config["g_diamond_blast_module_vf"]
    log:
        out = "log/linear_non_redundant_gene/gene_annotation/annotation_vf.out",
        err = "log/linear_non_redundant_gene/gene_annotation/annotation_vf.err"
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

rule liner_annotation_other_gene:
    input:
        f = "linear_non_redundant_gene/linear_gene_prodigal_protein_seq.faa"
    output:
        f = "linear_non_redundant_gene/annotation_BacMet2/annotation_BacMet2.txt"
    threads: config["threads"]
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    params:
        db = config["BacMet2_database"],
        outfmt = config["g_diamond_outfmt_BacMet2"],
        evalue = config["g_diamond_evalue_BacMet2"],
        id = config["g_diamond_id_BacMet2"],
        cover = config["g_diamond_cover_BacMet2"],
        md = config["g_diamond_blast_module_BacMet2"]
    log:
        out = "log/linear_non_redundant_gene/annotation_BacMet2.out",
        err = "log/linear_non_redundant_gene/annotation_BacMet2.err"
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
