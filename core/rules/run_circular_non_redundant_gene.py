#######################################################################################################
# !/usr/bin/env python3
# -*- coding:utf-8 -*-
# pipeline for plasmid and mobile
# 2020-0810
# wanli he
###########################################################################################################

# rule circular_genecalling_plasmid_gene:
#     input:
#         f = "circular/plasmid/{sample}_verify_plasmid_circular.fasta"
#     output:
#         f = temp("circular_non_redundant_gene/{sample}_nucl_gene.fa")
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

rule circular_genecalling_plasmid_gene:
    input:
        f = "circular_non_redundant_contig/circular_non_redundant_contigs"
    output:
        f = temp("circular_non_redundant_gene/circular_non_redundant_contigs_nucl_gene.fa")
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    params:
        p = config['g_prodigal_p']
    log:
        out = "log/linear_non_redundant_gene/circular_non_redundant_contigs_genecalling_genecalling.out",
        err = "log/linear_non_redundant_gene/circular_non_redundant_contigs_genecalling_genecalling.err"
    shell:
        "prodigal -p {params.p} " \
                "-i {input.f}_rep_seq.fasta " \
                "-d {output.f} " \
                "2>{log.err} >{log.out}"

rule rename_genecalling_file_cc:
    input:
        f = "circular_non_redundant_gene/circular_non_redundant_contigs_nucl_gene.fa",
    output:
        f = temp("circular_non_redundant_gene/circular_non_redundant_contigs_nucl_gene_rename.fa")
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

# rule circular_cutting_all_plasmid_gene:
#     input:
#         f = expand("circular_non_redundant_gene/{sample}_nucl_gene_rename.fa", sample=config["samples"])
#     output:
#         f =temp("circular_non_redundant_gene/all_plasmid_genes.fa")
#     shell:
#         "cat {input.f} >>{output.f}"


rule circular_cdhit_nucler_plasmid_gene:
    input:
        f1 = "circular_non_redundant_gene/circular_non_redundant_contigs_nucl_gene_rename.fa"
    output:
        f1 = "circular_non_redundant_gene/circular_non_redundant_genes.fa"
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        out = "log/circular_non_redundant_gene/cdhit_plasmid.out",
        err = "log/circular_non_redundant_gene/cdhit_plasmid.err"
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

rule circular_index_bam_plasmid_gene:
    input:
        f = "circular_non_redundant_gene/circular_non_redundant_genes.fa"
    output:
        f1= temp("circular_non_redundant_gene/circular_non_redundant_genes.fa.pac"),
        f2= temp("circular_non_redundant_gene/circular_non_redundant_genes.fa.amb"),
        f3= temp("circular_non_redundant_gene/circular_non_redundant_genes.fa.bwt"),
        f4= temp("circular_non_redundant_gene/circular_non_redundant_genes.fa.sa"),
        f5= temp("circular_non_redundant_gene/circular_non_redundant_genes.fa.ann")
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        out = "log/circular_index_bam_plasmid_gene/index.out",
        err = "log/circular_index_bam_plasmid_gene/index.err"
    shell:
        "bwa index {input.f};" \
        "rm -rf circular_non_redundant_gene/circular_non_redundant_genes.fa.clstr"


rule circular_get_bam_file_plasmid_gene:
    input:
        f1= "circular_non_redundant_gene/circular_non_redundant_genes.fa",
        f2= "qc_reads/{sample}_qc_1.fastq.gz",
        f3= "qc_reads/{sample}_qc_2.fastq.gz",
        f4= "circular_non_redundant_gene/circular_non_redundant_genes.fa.pac",
        f5= "circular_non_redundant_gene/circular_non_redundant_genes.fa.amb",
        f6= "circular_non_redundant_gene/circular_non_redundant_genes.fa.bwt",
        f7= "circular_non_redundant_gene/circular_non_redundant_genes.fa.sa",
        f8= "circular_non_redundant_gene/circular_non_redundant_genes.fa.ann"
    output:
        f=temp("circular_non_redundant_gene/{sample}_plasmid.bam")
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        out = "log/circular_get_bam_file_plasmid_gene/{sample}.out",
        err = "log/circular_get_bam_file_plasmid_gene/{sample}.err"
    shell:
        "bwa mem -M {input.f1} {input.f2} {input.f3} -t {threads}" \
        "| samtools view -Sb - >{output.f}"


rule circular_get_sort_file_plasmid_gene:
    input:
        f="circular_non_redundant_gene/{sample}_plasmid.bam"
    output:
        f= temp("circular_non_redundant_gene/{sample}_plasmid_sort.bam")
    threads: 10
    params:
        q= config["bwa_q"]
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    shell:
        "samtools sort -n -@ {threads} {input.f} -o {output.f}"

rule filter_msamtools_circular:
    input:
        f = "circular_non_redundant_gene/{sample}_plasmid_sort.bam"
    output:
        f =temp("circular_non_redundant_gene/{sample}_plasmid_sort_filter.bam")
    params:
        m = config["msamtools_path"],
        l = config["circular_msamtools_gene_l"],
        p = config["circular_msamtools_gene_p"],
        z = config["circular_msamtools_gene_z"]
    shell:
        "{params.m} filter -b -l {params.p} -p {params.p} -z {params.z} --besthit {input.f} >{output.f}"

rule circular_geneabundance:
    input:
        f = "circular_non_redundant_gene/{sample}_plasmid_sort_filter.bam"
    output:
        f = temp("circular_non_redundant_gene/{sample}.profile.txt")
    params:
        number = config["circular_gene_abundance_format"],
        m = config["msamtools_path"]
    run:
        id = os.path.basename(input.f)[:-len("_plasmid_sort_filter.bam")]
        os.system(params.m+" profile --multi=all  --label="+id+"  -o "+output.f+" "+input.f)


rule cutting_all_circular_profile_file:
    input:
        f = expand("circular_non_redundant_gene/{sample}.profile.txt",sample=config["samples"])
    output:
        f = "circular_non_redundant_gene/relative_gene_abundance/all_samples_gene_relative_abundance.txt"
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



#gene functional annotation
rule circular_gene_genecalling_gene:     #non-redundant-gene-set
    input:
        f = "circular_non_redundant_gene/circular_non_redundant_genes.fa"
    output:
        f = temp("circular_non_redundant_gene/gene_prodigal_protein_seq__.faa")
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    params:
        p = config["g_prodigal_p"]
    log:
        out ="log/circular_gene_annotation/plasmid_gene_genecalling.out",
        err ="log/circular_gene_annotation/plasmid_gene_genecalling.err"
    shell:
        "prodigal -p {params.p} " \
                "-i {input.f} " \
                "-a {output.f} " \
                "2>{log.err} >{log.out}"

rule rename_genecalling_file_cc_aa:
    input:
        f = "circular_non_redundant_gene/gene_prodigal_protein_seq__.faa",
    output:
        f = temp("circular_non_redundant_gene/gene_prodigal_protein_seq.faa")
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


rule circular_functional_annotation_genes:
    input:
        f = "circular_non_redundant_gene/circular_non_redundant_genes.fa"
    output:
        f = directory("circular_non_redundant_gene/functional_annotation")
    threads: config['functional_threads']
    conda:
        "%s/emapper.yaml" % CONDAENV
    log:
        out = "log/circular_gene_annotation/functional_annotation.out",
        err = "log/circular_gene_annotation/functional_annotation.err"
    params:
        m = config['g_functional_annotation_m'],
        d = E_DB
    script:
        "../scripts/emapper.py"



#other gene annotation
rule circular_annotation_ARGs_gene:
    input:
        f = "circular_non_redundant_gene/circular_non_redundant_genes.fa"
    output:
        f = directory("circular_non_redundant_gene/annotation_ARGs")
    threads: config["threads"]
    conda:
        "%s/rgi.yaml" % CONDAENV
    log:
        out = "log/circular_gene_annotation/annotation_ARGs.out",
        err = "log/circular_gene_annotation/annotation_ARGs.err"
    params:
        d = config['g_rgi_d'],
        it = config['g_rgi_input_type'],
        a = config['g_rgi_a'],
        mode = config["g_rgi_mode"],
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


rule circular_annotation_vf_gene:
    input:
        f = "circular_non_redundant_gene/gene_prodigal_protein_seq.faa"
    output:
        f = "circular_non_redundant_gene/annotation_vf/annotation_vf.txt"
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
        out = "log/circular_non_redundant_gene/gene_annotation/annotation_vf.out",
        err = "log/circular_non_redundant_gene/gene_annotation/annotation_vf.err"
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

rule circular_annotation_other_gene:
    input:
        f = "circular_non_redundant_gene/gene_prodigal_protein_seq.faa"
    output:
        f = "circular_non_redundant_gene/annotation_BacMet2/annotation_BacMet2.txt"
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
        out = "log/circular_non_redundant_gene/annotation_BacMet2.out",
        err = "log/circular_non_redundant_gene/annotation_BacMet2.err"
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

# rule circular_final_summary:
#     input:
#         f1 = "circular_non_readundant_gene/relative_gene_abundance/all_samples_gene_relative_abundance.txt",
#         f2 = "circular_non_readundant_gene/annotation_ARGs",
#         f3 = "circular_non_readundant_gene/annotation_vf/annotation_vf.txt",
#         f4 = "circular_non_readundant_gene/annotation_BacMet2/annotation_BacMet2.txt"
#     output:
#         f = "summary_file/circular_gene_suammary_report.txt"
#     run:
#         (dict_ab,d_arg,d_vf,d_bac)=({},{},{},{})
#         with open(input.f1,"r") as infile1:
#             with open(input.f2,"r") as infile2:
#                 with open(input.f3,"r") as infile3:
#                     with open(input.f4,"r") as infile4:
#                         header_ab = infile1.readline()
#                         for line in infile1:
#                             lst=line.strip("\t",1)
#                             dict_ab[lst[0]] = lst[1]
#                             d_arg[lst[0]] = "-"
#                             d_vf[lst[0]] = "-"
#                             d_bac[lst[0]] = "-"
#
#                         header_args=infile2.readline()
#                         for line in infile2:
#                             lst2=line.strip().split()
#                             if lst2[0] in d_arg.keys():
#                                 d_arg[lst2[0]] = lst2[16]
#                             else:
#                                 pass
#
#                         header_vf = infile3.readline()
#                         for line in infile3:
#                             lst3=line.split().split()
#                             if lst3[0] in d_vf.keys():
#                                 d_vf[lst3[0]] = lst3[]   ######
#                             else:
#                                 pass
#
#                         header_bac = infile4.readline()
#                         for line in infile4:
#                             lst4 = line.strip().split()
#                             if lst4[0] in d_bac.keys():
#                                 d_bac[lst4[0]] = lst4[]  ######
#                             else:
#                                 pass






























