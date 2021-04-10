#######################################################################################################
# !/usr/bin/env python3
# -*- coding:utf-8 -*-
# pipeline for plasmid and mobile
# 2019-03-06
# wanli he
###########################################################################################################


# rule linear_genecalling_plasmid_gene:
#     input:
#         f = "linear_plasmid_genome/{sample}_predict_plasmid.fa"
#     output:
#         f = temp("linear_non_redundant_gene/{sample}_nucl_gene.fa")
#     threads: config['threads']
#     conda:
#         "%s/non_redundant.yaml" % CONDAENV
#     params:
#         p = config['g_prodigal_p']
#     log:
#         out = "log/linear_non_redundant_gene/{sample}_genecalling_genecalling.out",
#         err = "log/linear_non_redundant_gene/{sample}_genecalling_genecalling.err"
#     script:
#         "../scripts/gene_calling.py"

rule linear_genecalling_plasmid_gene:
    input:
        f = "linear_non_redundant_contig/linear_non_redundant_contigs"
    output:
        f = temp("linear_non_redundant_gene/linear_non_redundant_contigs_nucl_gene.fa")
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    params:
        p = config['g_prodigal_p']
    log:
        out = "log/linear_non_redundant_gene/linear_non_redundant_contigs_genecalling_genecalling.out",
        err = "log/linear_non_redundant_gene/linear_non_redundant_contigs_genecalling_genecalling.err"
    shell:
        "prodigal -p {params.p} " \
                "-i {input.f}_rep_seq.fasta " \
                "-d {output.f} " \
                "2>{log.err} >{log.out}"


rule rename_genecalling_file_ff:
    input:
        f = "linear_non_redundant_gene/linear_non_redundant_contigs_nucl_gene.fa",
    output:
        f = temp("linear_non_redundant_gene/linear_non_redundant_contigs_nucl_gene_rename.fa")
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

# rule linear_cutting_all_plasmid_gene:
#     input:
#         f = expand("linear_non_redundant_gene/{sample}_nucl_gene_rename.fa", sample=config["samples"])
#     output:
#         f =temp("linear_non_redundant_gene/all_plasmid_genes.fa")
#     shell:
#         "cat {input.f} >>{output.f}"


rule linear_cdhit_nucler_plasmid_gene:
    input:
        f1 = "linear_non_redundant_gene/linear_non_redundant_contigs_nucl_gene_rename.fa"
    output:
        f1 = "linear_non_redundant_gene/linear_non_redundant_genes.fa"
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        out = "log/linear_non_redundant_gene/cdhit_plasmid.out",
        err = "log/linear_non_redundant_gene/cdhit_plasmid.err"
    params:
        cd = config['cdhit-est_path'],
        c = config['g_cd-hit_c'],
        aS = config['g_cd-hit_aS'],
        M = config['g_cd-hit_M']
    shell:
        "{params.cd} -i {input.f1} " \
                    "-c {params.c} " \
                    "-M {params.M} " \
                    "-aS {params.aS} " \
                    "-o {output.f1} " \
                    "-T {threads}" \
                    " 2>{log.err} >{log.out}"

rule linear_index_bam_plasmid_gene:
    input:
        f = "linear_non_redundant_gene/linear_non_redundant_genes.fa"
    output:
        f1= temp("linear_non_redundant_gene/linear_non_redundant_genes.fa.pac"),
        f2= temp("linear_non_redundant_gene/linear_non_redundant_genes.fa.amb"),
        f3= temp("linear_non_redundant_gene/linear_non_redundant_genes.fa.bwt"),
        f4= temp("linear_non_redundant_gene/linear_non_redundant_genes.fa.sa"),
        f5= temp("linear_non_redundant_gene/linear_non_redundant_genes.fa.ann")
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        out = "log/linear_index_bam_plasmid_gene/index.out",
        err = "log/linear_index_bam_plasmid_gene/index.err"
    shell:
        "bwa index {input.f};"
        "rm -rf linear_non_redundant_gene/linear_non_redundant_genes.fa.clstr"


rule linear_get_bam_file_plasmid_gene:
    input:
        f1= "linear_non_redundant_gene/linear_non_redundant_genes.fa",
        f2= "qc_reads/{sample}_qc_1.fastq.gz",
        f3= "qc_reads/{sample}_qc_2.fastq.gz",
        f4= "linear_non_redundant_gene/linear_non_redundant_genes.fa.pac",
        f5= "linear_non_redundant_gene/linear_non_redundant_genes.fa.amb",
        f6= "linear_non_redundant_gene/linear_non_redundant_genes.fa.bwt",
        f7= "linear_non_redundant_gene/linear_non_redundant_genes.fa.sa",
        f8= "linear_non_redundant_gene/linear_non_redundant_genes.fa.ann"
    output:
        f=temp("linear_non_redundant_gene/{sample}_plasmid.bam")
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        out = "log/circular_get_bam_file_plasmid_gene/{sample}.out",
        err = "log/circular_get_bam_file_plasmid_gene/{sample}.err"
    shell:
        "bwa mem -M {input.f1} {input.f2} {input.f3} -t {threads}" \
        "| samtools view -Sb - >{output.f}"


rule linear_get_sort_file_plasmid_gene:
    input:
        f="linear_non_redundant_gene/{sample}_plasmid.bam"
    output:
        f=temp("linear_non_redundant_gene/{sample}_plasmid_sort.bam")
    threads: config['threads']
    params:
        q = config["bwa_q"]
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    shell:
        "samtools sort -n -@ {threads} {input.f} -o {output.f}"

rule filter_msamtools_linear:
    input:
        f = "linear_non_redundant_gene/{sample}_plasmid_sort.bam"
    output:
        f =temp("linear_non_redundant_gene/{sample}_plasmid_sort_filter.bam")
    params:
        m = config["msamtools_path"],
        l = config["linear_msamtools_gene_l"],
        p = config["linear_msamtools_gene_p"],
        z = config["linear_msamtools_gene_z"]
    shell:
        "{params.m} filter -b -l {params.p} -p {params.p} -z {params.z} --besthit {input.f} >{output.f}"

rule linear_geneabundance:
    input:
        f = "linear_non_redundant_gene/{sample}_plasmid_sort_filter.bam"
    output:
        f = temp("linear_non_redundant_gene/{sample}.profile.txt")
    params:
        # number = config["circular_gene_abundance_format"],
        m = config["msamtools_path"]
    run:
        id = os.path.basename(input.f)[:-len("_plasmid_sort_filter.bam")]
        os.system(params.m+" profile --multi=all  --label="+id+"  -o "+output.f+" "+input.f)

rule cutting_all_linear_profile_file:
    input:
        f = expand("linear_non_redundant_gene/{sample}.profile.txt",sample=config["samples"])
    output:
        f = "linear_non_redundant_gene/relative_gene_abundance/all_samples_gene_relative_abundance.txt"
    run:
        if os.path.exists(output.f):
            pass
        else:
            dict = defaultdict(str)
            header = "gene_ID"

            for file in input.f:
                sample_ID = os.path.basename(file)[:-len(".profile.txt")]
                head = "\t" + sample_ID
                header += head
                with open(file, "r") as infile:
                    infile.readline()
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

# rule linear_reads_mapping_gene:
#     input:
#         f="linear_non_readundant_gene/{sample}_plasmid_sort.bam"
#     output:
#         f=temp("linear_non_readundant_gene/mapping/{sample}_bp_mapping.txt")
#     threads: config['threads']
#     conda:
#         "%s/bedtools.yaml" % CONDAENV
#     shell:
#         "genomeCoverageBed -ibam {input.f} -d > {output.f}"

# rule filter_dectect_linear_gene:
#     input:
#         f = "linear_non_readundant_gene/mapping/{sample}_bp_mapping.txt",
#         f2 = "linear_non_readundant_gene/linear_non_redundant_genes.fa"
#     output:
#         f = temp("linear_non_readundant_gene/{sample}_gene_coverage.txt")
#     run:
#         dict_len={}
#         with open(input.f2,"r") as infile2:
#             for line in infile2:
#                 if line.startswith(">"):
#                     k = line.strip().split()[0][1:]
#                     dict_len[k] = ""
#                 else:
#                     dict_len[k] += line.strip()
#
#         dict_c = defaultdict(float)
#         dict_det = defaultdict(float)
#         with open(input.f,"r") as infile:
#             for line in infile:
#                 lst = line.strip().split()
#                 dict_c[lst[0]] += float(lst[2])
#
#                 if float(lst[2]) != 0:
#                     dict_det[lst[0]] += 1
#                 else:
#                     dict_det[lst[0]] += 0
#
#         with open(output.f,"w") as outfile:
#             for k in dict_det.keys():
#                 t = dict_det[k]/float(len(dict_len[k]))
#                 if t > float(config["gene_detection"]):
#                     st = "{}\t{}\n".format(k,dict_c[k])
#                     outfile.write(st)
#                 else:
#                     st = "{}\t{}\n".format(k,0)
#                     outfile.write(st)
#
#
# rule linear_relative_gene_abundance_gene:
#     input:
#         f = "linear_non_readundant_gene/{sample}_gene_coverage.txt"
#     output:
#         f = temp("linear_non_readundant_gene/{sample}_relative_abundance.txt")
#     run:
#         a = 0
#         with open(input.f,"r") as infile:
#             for line in infile:
#                 lst = line.strip().split()
#                 a += float(lst[1])
#
#         with open(input.f,"r") as infile:
#             with open(output.f,"w") as outfile:
#                 for line in infile:
#                     lst = line.strip().split()
#                     if a != 0:
#                         s= float(lst[1])/float(a)
#                         st = lst[0] +"\t" + str(s) +"\n"
#                         outfile.write(st)
#                     else:
#                         st = lst[0] +"\t" +"0\n"
#                         outfile.write(st)
#
# rule linear_paste_relative_gene_abundance_gene:
#     input:
#         f = expand("linear_non_readundant_gene/{sample}_relative_abundance.txt",sample=config["samples"])
#     output:
#         f = "linear_non_readundant_gene/relative_gene_abundance/all_samples_gene_relative_abundance.txt"
#     run:
#         if os.path.exists(output.f):
#             pass
#         else:
#             dict = defaultdict(str)
#             header = "gene_ID"
#
#             for file in input.f:
#                 sample_ID = os.path.basename(file)[:-len("_relative_abundance.txt")]
#                 head = "\t" + sample_ID
#                 header += head
#                 with open(file, "r") as infile:
#                     for line in infile:
#                         lst = line.strip().split()
#                         lt = "\t" + lst[1]
#                         dict[lst[0]] += lt
#
#             with open(output.f, "w") as outfile:
#                 header += "\n"
#                 outfile.write(header)
#                 for k in dict.keys():
#                     st = "{}\t{}\n".format(k, dict[k])
#                     outfile.write(st)
#


#gene functional annotation
rule linear_gene_genecalling_gene:     #non-redundant-gene-set
    input:
        f = "linear_non_redundant_gene/linear_non_redundant_genes.fa"
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
                "-i {input.f} " \
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
        f = "linear_non_redundant_gene/linear_non_redundant_genes.fa"
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
        f = "linear_non_redundant_gene/linear_non_redundant_genes.fa"
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
        "rgi main -i {input.f} " \
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
