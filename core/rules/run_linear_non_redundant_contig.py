#######################################################################################################
# !/usr/bin/env python3
# -*- coding:utf-8 -*-
# pipeline for plasmid and mobile
# 2019-03-06
# wanli he
###########################################################################################################



# contig abundance
rule cutting_linear_plasmid:
    input:
        f = expand("linear_plasmid_genome/{sample}_predict_plasmid.fa", sample=config["samples"])
    output:
        f =temp("linear_non_redundant_contig/all_plasmid_contigs.fa")
    shell:
        "cat {input.f} >>{output.f}"


# rule cdhit_nucler_linear_plasmid:
#     input:
#         f1 = "linear_non_redundant_contig/all_plasmid_contigs.fa"
#     output:
#         f1 = "linear_non_redundant_contig/linear_non_redundant_contigs.fa"
#     threads: config['threads']
#     conda:
#         "%s/non_redundant.yaml" % CONDAENV
#     log:
#         out = "log/linear_non_redundant_contig/cdhit_plasmid.out",
#         err = "log/linear_non_redundant_contig/cdhit_plasmid.err"
#     params:
#         cd = config['cdhit-est_path'],
#         c = config['c_cd-hit_c'],
#         aS = config['c_cd-hit_aS'],
#         M = config['c_cd-hit_M']
#     shell:
#         "{params.cd} -i {input.f1} " \
#                     "-c {params.c} " \
#                     "-M {params.M} " \
#                     "-aS {params.aS} " \
#                     "-o {output.f1} " \
#                     "-T {threads}" \
#                     " 2>{log.err} >{log.out}"

rule non_redundant_nucler_circular_plasmid:
    input:
        f = "linear_non_redundant_contig/all_plasmid_contigs.fa"
    output:
        f = "linear_non_redundant_contig/linear_non_redundant_contigs"
    threads: config['threads']
    conda:
        "%s/mmseqs2.yaml" % CONDAENV
    log:
        out = "log/linear_non_redundant_contig/non_redundant_plasmid.out",
        err = "log/linear_non_redundant_contig/non_redundant_plasmid.err"
    params:
        id = config["linear_contig_mmseqs2_id"],
        c =  config["linear_contig_mmseqs2_c"],
        mode =  config["linear_contig_mmseqs2_mode"]
    shell:
        "mmseqs easy-cluster {input.f} {output.f} tmp --min-seq-id {params.id} -c {params.c} " \
        "--cov-mode {params.mode} 2>{log.err} >{log.out};" \
        "touch {output.f}"


# rule linear_plasmidverify:
#     input:
#         f1="linear_non_redundant_contig/linear_non_redundant_contigs"
#     output:
#         f2=directory("linear_non_redundant_contig/Non_redundant_contig_plasmidverify")
#     threads: config['linear_plasmidverify_threads']
#     conda:
#         "%s/plasmidverify.yaml" % CONDAENV
#     params:
#         hmm=config['Pfam_database']
#     log:
#         out = "log/linear_plasmidverify/plasmidverify_liner.out",
#         err = "log/linear_plasmidverify/plasmidverify_liner.err"
#     shell:
#         "{config[plasmidverify_path]} -f {input.f1}_rep_seq.fasta " \
#                              "--hmm {params.hmm} " \
#                              "-t {threads} " \
#                              "-o {output.f2}" \
#                              " 2>{log.err} >{log.out}"


rule linear_index_bam_plasmid:
    input:
        f = "linear_non_redundant_contig/linear_non_redundant_contigs"
    output:
        f1= temp("linear_non_redundant_contig/linear_non_redundant_contigs_rep_seq.fasta.pac"),
        f2= temp("linear_non_redundant_contig/linear_non_redundant_contigs_rep_seq.fasta.amb"),
        f3= temp("linear_non_redundant_contig/linear_non_redundant_contigs_rep_seq.fasta.bwt"),
        f4= temp("linear_non_redundant_contig/linear_non_redundant_contigs_rep_seq.fasta.sa"),
        f5= temp("linear_non_redundant_contig/linear_non_redundant_contigs_rep_seq.fasta.ann")
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        err = "log/linear_index_bam_plasmid/index.err",
        out = "log/linear_index_bam_plasmid/index.out"
    shell:
        "bwa index {input.f}_rep_seq.fasta ;"
        "rm -rf linear_non_redundant_contig/linear_non_redundant_contigs_rep_seq.fasta.clstr"



rule linear_get_bam_file_plasmid:
    input:
        f1= "linear_non_redundant_contig/linear_non_redundant_contigs",
        f2= "qc_reads/{sample}_qc_1.fastq.gz",
        f3= "qc_reads/{sample}_qc_2.fastq.gz",
        f4= "linear_non_redundant_contig/linear_non_redundant_contigs_rep_seq.fasta.pac",
        f5= "linear_non_redundant_contig/linear_non_redundant_contigs_rep_seq.fasta.amb",
        f6= "linear_non_redundant_contig/linear_non_redundant_contigs_rep_seq.fasta.bwt",
        f7= "linear_non_redundant_contig/linear_non_redundant_contigs_rep_seq.fasta.sa",
        f8= "linear_non_redundant_contig/linear_non_redundant_contigs_rep_seq.fasta.ann"
    output:
        f= temp("linear_non_redundant_contig/{sample}_plasmids.bam")
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        err = "log/linear_get_bam_file_plasmid/{sample}.err",
        out = "log/linear_get_bam_file_plasmid/{sample}.out"
    shell:
        "bwa mem -M {input.f1}_rep_seq.fasta {input.f2} {input.f3} -t {threads} | samtools view -Sb - >{output.f}"


rule linear_get_sort_file_plasmid:
    input:
        f="linear_non_redundant_contig/{sample}_plasmids.bam"
    output:
        f=temp("linear_non_redundant_contig/{sample}_plasmids_sort.bam")
    threads: config['threads']
    params:
        q= config["bwa_q"]
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    shell:
        "samtools view -h {input.f}|" \
        "samtools view -Su -q{params.q} -|samtools sort -O BAM -o {output.f} -"


rule linear_reads_mapping:
    input:
        f="linear_non_redundant_contig/{sample}_plasmids_sort.bam"
    output:
        f=temp("linear_non_redundant_contig/{sample}_bp_mapping.txt")
    threads: config['threads']
    conda:
        "%s/bedtools.yaml" % CONDAENV
    shell:
        "genomeCoverageBed -ibam {input.f} -d > {output.f}"


rule filter_dectect_linear_contig:
    input:
        f = "linear_non_redundant_contig/{sample}_bp_mapping.txt",
        f2 = "linear_non_redundant_contig/linear_non_redundant_contigs"
    output:
        f = temp("linear_non_redundant_contig/{sample}_coverage.txt")
    run:
        dict_len={}
        with open(input.f2+'_rep_seq.fasta',"r") as infile2:
            for line in infile2:
                if line.startswith(">"):
                    k = line.strip().split()[0][1:]
                    dict_len[k] = ""
                else:
                    dict_len[k] += line.strip()

        dict_c = defaultdict(float)
        dict_det = defaultdict(float)
        with open(input.f,"r") as infile:
            for line in infile:
                lst = line.strip().split()
                dict_c[lst[0]] += float(lst[2])

                if float(lst[2]) != 0:
                    dict_det[lst[0]] += 1
                else:
                    dict_det[lst[0]] += 0

        with open(output.f,"w") as outfile:
            for k in dict_det.keys():
                t = dict_det[k]/float(len(dict_len[k]))
                if t > float(config["contig_detection"]):
                    st = "{}\t{}\n".format(k,dict_c[k])
                    outfile.write(st)
                else:
                    st = "{}\t{}\n".format(k,0)
                    outfile.write(st)



rule linear_relative_contig_abundance:
    input:
        f = "linear_non_redundant_contig/{sample}_coverage.txt"
    output:
        f = temp("linear_non_redundant_contig/{sample}_relative_abundance.txt")
    threads: config['threads']
    run:
        # data = pd.read_csv(input.f,sep="\t")
        # a = data.iloc[:,1].sum()
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

rule linear_paste_relative_contig_abundance:
    input:
        f = expand("linear_non_redundant_contig/{sample}_relative_abundance.txt",sample=config["samples"])
    output:
        f = "linear_non_redundant_contig/relative_contig_abundance/all_samples_contig_relative_abundance.txt"
    threads: config['threads']
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
                        lt = lst[1]
                        dict[lst[0]] += lt

            with open(output.f, "w") as outfile:
                header += "\n"
                outfile.write(header)
                for k in dict.keys():
                    st = "{}\t{}\n".format(k, dict[k])
                    outfile.write(st)

# classify
# if config["linearized_tool"] != "platon":
#     rule linear_classify_platon:
#         input:
#             f = "linear_non_readundant_contig/linear_non_redundant_contigs.fa"
#         output:
#             f ="linear_non_readundant_contig/clasify_platon"
#         threads: config["classify_platon_threads"]
#         conda:
#             "%s/linearized-platon.yaml" % CONDAENV
#         params:
#             mode = config["classify_platon_mode"],
#             db = config["platondb"]
#         log:
#             out = "log/classify-linear-contig/platon_linear.out",
#             err = "log/classify-linear-contig/platon_linear.err"
#         shell:
#             "platon {input.f} " \
#                     "--db {params.db} " \
#                     "--mode {params.mode} " \
#                     "-c -o {output.f} " \
#                     "-t {threads} " \
#                     "-v 2>{log.err} >{log.out}"


# rule linear_split_contigs:
#     input:
#         f1="linear_non_readundant_contig/linear_non_redundant_contigs.fa"
#     output:
#         f2=temp(directory("linear_plasmid_contig_split"))
#     threads: config["threads"]
#     run:
#         if not os.path.exists(output.f2):
#             os.mkdir(output.f2)
#         with open(input.f1,mode='r') as infile:
#             dir = defaultdict(str)
#             for line in infile:
#                 if line.startswith(">"):
#                     k = line
#                     dir[k] = ""
#                 else:
#                     dir[k] += line
#
#             for k, v in dir.items():
#                 res = k[1:].split()[0] + "_split.fa"
#                 a = os.path.join(output.f2, res)
#                 with open(a,mode="w") as outfile:
#                     outfile.write(k)
#                     outfile.write(v)
#
#
# rule linear_run_mob_typer:
#     input:
#         f1="linear_plasmid_contig_split"
#     output:
#         f2=temp(directory("linear_non_readundant_contig/all_result_mob_typer"))
#     threads: config["threads_mob_typer"]
#     conda:
#         "%s/mobtyper.yaml" % CONDAENV
#     log:
#         err="log/linear_run_classify_plasmid/run_mob_typer.err",
#         out="log/linear_run_classify_plasmid/run_mob_typer.out"
#     script:
#         "../scripts/run_mob_typer.py"
#
#
#
# rule linear_clean_mob_typer_file:
#     input:
#         f1="linear_non_readundant_contig/all_result_mob_typer"
#     output:
#         f2="linear_non_readundant_contig/plasmid_classify/plasmid_classify.txt"
#     threads: config["threads"]
#     run:
#         path=input.f1
#         files=glob.glob(path+"/*_mob_typer_result")
#         lst=[]
#         header=""
#         for file in files:
#             son_path = file+"/*"
#             for son_file in glob.glob(son_path):
#                 with open(son_file,mode="r") as infile:
#                     header = infile.readline()
#                     st=infile.readline()
#                     lst.append(st)
#
#         with open(output.f2,mode="w") as outfile:
#             outfile.write(header)
#             for ele in lst:
#                 outfile.write(ele)
#

#co-exist
rule linear_makefaa:
    input:
        f="linear_non_redundant_contig/linear_non_redundant_contigs"
    output:
        f1 = temp("linear_non_redundant_contig/all_genecalling_protein__.faa"),
        f2 = temp("linear_non_redundant_contig/all_genecalling_nucl__.fa")
    threads: config["threads"]
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        out = "log/linear_co-exist_circular/all_genecalling.out",
        err = "log/linear_co-exist_circular/all_genecalling.err"
    params:
        p = config["c_prodigal_p"]
    shell:
        "prodigal -p {params.p} -i {input.f}_rep_seq.fasta -a {output.f1} -d {output.f2} 2>{log.err} >{log.out}"

rule rename_genecalling_file_vvv:
    input:
        f = "linear_non_redundant_contig/all_genecalling_protein__.faa",
        f2 = "linear_non_redundant_contig/all_genecalling_nucl__.fa"
    output:
        f = temp("linear_non_redundant_contig/all_genecalling_protein.faa"),
        f2 = temp("linear_non_redundant_contig/all_genecalling_nucl.fa")
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

        with open(input.f2,"r") as infile:
            with open({output.f2},"w") as outfile:
                for line in infile:
                    if line.startswith(">"):
                        lst = line.strip().split("\t",1)
                        la = lst[0].rsplit("_",1)[0]
                        st = la+"\n"
                        outfile.write(st)
                    else:
                        outfile.write(line)

rule linear_plasmid_with_MGEs:
    input:
        f = "linear_non_redundant_contig/all_genecalling_protein.faa"
    output:
        f = "linear_non_redundant_contig/co-exist/all_plasmids_with_MEGs.txt"
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


rule linear_contig_annotation_ARGs:
    input:
        f = "linear_non_redundant_contig/linear_non_redundant_contigs"
    output:
        f = directory("linear_non_redundant_contig/co-exist/annotation_ARGs")
    threads: config["threads"]
    conda:
        "%s/rgi.yaml" % CONDAENV
    log:
        out = "log/co-exist/linear_ARGs.out",
        err = "log/co-exist/linear_ARGs.err"
    params:
        d = config['c_rgi_d'],
        it = config['c_rgi_input_type'],
        a = config['c_rgi_a'],
        mode = config["c_rgi_mode"]
    shell:
        "mkdir -p {output.f};" \
        "rgi main -i {input.f}_rep_seq.fasta" \
                "--output_file {output.f} " \
                "-d {params.d} " \
                "--input_type {params.it} " \
                "-n {threads} " \
                "-a {params.a} " \
                "{params.mode} --clean " \
                "2>{log.err} >{log.out};" \
        "mv {output.f}.json {output.f}.txt {output}"



rule linear_contig_annotation_vf:
    input:
        f = "linear_non_redundant_contig/all_genecalling_protein.faa"
    output:
        f = "linear_non_redundant_contig/co-exist/annotation_vf/annotation_vf.txt"
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
        out = "log/co-exist/linear_vf.out",
        err = "log/co-exist/linear_vf.err"
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

rule linear_contig_annotation_other:
    input:
        f = "linear_non_redundant_contig/all_genecalling_protein.faa"
    output:
        f = "linear_non_redundant_contig/co-exist/annotation_BacMet2/annotation_BacMet2.txt"
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
        out = "log/co-exist/linear_BacMet2.out",
        err = "log/co-exist/linear_BacMet2.err"
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

# rule linear_contig_annotation_plsdb:
#     input:
#         f = "linear_non_readundant_contig/linear_non_redundant_contigs.fa"
#     output:
#         f = "linear_non_readundant_contig/co-exist/annotation_plsdb/annotation_plsdb.txt"
#     threads: config["threads"]
#     conda:
#         "%s/non_redundant.yaml" % CONDAENV
#     params:
#         db = config["plsdb_database"],
#         outfmt = config["c_diamond_outfmt_plsdb"],
#         evalue = config["c_diamond_evalue_plsdb"],
#         id = config["c_diamond_id_plsdb"],
#         cover = config["c_diamond_cover_plsdb"],
#         md = config["c_diamond_blast_module_plsdb"]
#     log:
#         out = "log/co-exist/linear_plsdb.out",
#         err = "log/co-exist/linear_plsdb.err"
#     shell:
#         "diamond {params.md} " \
#                 "--db {params.db} " \
#                 "--query {input.f} " \
#                 "--out {output.f} " \
#                 "--outfmt {params.outfmt} " \
#                 "--evalue {params.evalue} " \
#                 "--id {params.id} " \
#                 "--query-cover {params.cover} " \
#                 "--threads {threads} " \
#                 "2>{log.err} >{log.out}"

#gene synteny
rule linear_gene_synteny:
    input:
        f = "linear_non_redundant_contig/all_genecalling_nucl.fa"
    output:
        f = "linear_non_redundant_contig/gene_synteny.txt"
    run:
        dict = defaultdict(str)
        with open(input.f,"r") as infile:
            for line in infile:
                if line.startswith(">"):
                    lst = line.strip().split()
                    # if lst[8].split(";")[1] == "partial=00":
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