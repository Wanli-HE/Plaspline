#######################################################################################################
# !/usr/bin/env python3
# -*- coding:utf-8 -*-
# pipeline for plasmid and mobile
# 2021-1-21
# wanli he
###########################################################################################################

ruleorder:

rule cutting_circular_plasmid:
    input:
        f = expand("circular/plasmid/{sample}_verify_plasmid_circular.fasta", sample=config["samples"])
    output:
        f =temp("circular_non_redundant_contig/all_plasmid_contigs.fa")
    shell:
        "cat {input.f} >>{output.f}"

rule non_redundant_nucler_circular_plasmid__:
    input:
        f = "circular_non_redundant_contig/all_plasmid_contigs.fa"
    output:
        f = "circular_non_redundant_contig/circular_non_redundant_contigs"
    threads: config['threads']
    conda:
        "%s/mmseqs2.yaml" % CONDAENV
    log:
        out = "log/circular_non_redundant_contig/non_redundant_plasmid.out",
        err = "log/circular_non_redundant_contig/non_redundant_plasmid.err"
    params:
        id = config["circular_contig_mmseqs2_id"],
        c =  config["circular_contig_mmseqs2_c"],
        mode =  config["circular_contig_mmseqs2_mode"]
    shell:
        "mmseqs easy-cluster {input.f} {output.f} tmp --min-seq-id {params.id} -c {params.c} --cov-mode {params.mode} " \
        "2>{log.err} >{log.out};touch {output.f}"

rule circular_index_bam_plasmid:
    input:
        f = "circular_non_redundant_contig/circular_non_redundant_contigs"
    output:
        f1= temp("circular_non_redundant_contig/circular_non_redundant_contigs_rep_seq.fasta.pac"),
        f2= temp("circular_non_redundant_contig/circular_non_redundant_contigs_rep_seq.fasta.amb"),
        f3= temp("circular_non_redundant_contig/circular_non_redundant_contigs_rep_seq.fasta.bwt"),
        f4= temp("circular_non_redundant_contig/circular_non_redundant_contigs_rep_seq.fasta.sa"),
        f5= temp("circular_non_redundant_contig/circular_non_redundant_contigs_rep_seq.fasta.ann")
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        err = "log/circular_index_bam_plasmid/index.err",
        out = "log/circular_index_bam_plasmid/index.out"
    shell:
        "bwa index {input.f}_rep_seq.fasta 2>{log.err} >{log.out};"
        "rm -rf circular_non_redundant_contig/circular_non_redundant_contigs.fa.clstr"


rule circular_get_bam_file_plasmid:
    input:
        f1= "circular_non_redundant_contig/circular_non_redundant_contigs",
        f2= "qc_reads/{sample}_qc_1.fastq.gz",
        f3= "qc_reads/{sample}_qc_2.fastq.gz",
        f4= "circular_non_redundant_contig/circular_non_redundant_contigs_rep_seq.fasta.pac",
        f5= "circular_non_redundant_contig/circular_non_redundant_contigs_rep_seq.fasta.amb",
        f6= "circular_non_redundant_contig/circular_non_redundant_contigs_rep_seq.fasta.bwt",
        f7= "circular_non_redundant_contig/circular_non_redundant_contigs_rep_seq.fasta.sa",
        f8= "circular_non_redundant_contig/circular_non_redundant_contigs_rep_seq.fasta.ann"
    output:
        f= temp("circular_non_redundant_contig/{sample}_plasmids.bam")
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        err = "log/circular_get_bam_file_plasmid/{sample}.err",
        out = "log/circular_get_bam_file_plasmid/{sample}.out"
    shell:
        "bwa mem -M {input.f1}_rep_seq.fasta {input.f2} {input.f3} -t {threads} 2>{log.err}| samtools view -Sb - >{output.f}"


rule circular_get_sort_file_plasmid:
    input:
        f="circular_non_redundant_contig/{sample}_plasmids.bam"
    output:
        f=temp("circular_non_redundant_contig/{sample}_plasmids_sort.bam")
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        err = "log/circular_get_sort_file_plasmid/{sample}.err",
        out = "log/circular_get_sort_file_plasmid/{sample}.out"
    shell:
        "samtools sort -n -@ {threads} {input.f} -o {output.f} 2>{log.err} >{log.out}"

rule circular_reads_mapping:
    input:
        f="circular_non_redundant_contig/{sample}_plasmids_sort.bam"
    output:
        f=temp("circular_non_redundant_contig/{sample}_plasmid_sort_filter.bam")
    threads: config['threads']
    conda:
        "%s/msamtools.yaml" % CONDAENV
    params:
        l = config["circular_msamtools_l"],
        p = config["circular_msamtools_p"],
        z = config["circular_msamtools_z"]
    shell:
        "msamtools filter -b -l {params.l} -p {params.p} -z {params.z} --besthit {input.f} >{output.f}"

rule circular_gene_abundance_msa:
    input:
        f = "circular_non_redundant_contig/{sample}_plasmid_sort_filter.bam"
    output:
        f = temp("circular_non_redundant_contig/{sample}_profile_rb.txt")
    params:
        ff = config["circular_msamtools_contig_abundance_unit"],
    conda:
        "%s/msamtools.yaml" % CONDAENV
    shell:
        "msamtools profile --multi=all --unit={params.ff}  " \
        "--label={input.f} -o {output.f} {input.f}"


rule circular_coverage:
    input:
        f = "circular_non_redundant_contig/{sample}_plasmid_sort_filter.bam"
    output:
        f = temp("circular_non_redundant_contig/{sample}_coverage_info.txt.gz")
    conda:
        "%s/msamtools.yaml" % CONDAENV
    shell:
        "msamtools coverage -z --summary -o {output.f} {input.f}"


rule existing_circular_contigs:
    input:
        f = "circular_non_redundant_contig/{sample}_coverage_info.txt.gz",
        f2 = "circular_non_redundant_contig/{sample}_profile_rb.txt"
    output:
        f = "circular_non_redundant_contig/contig_abundance/per_sample/{sample}_existing_contigs_ab_per_sample.txt"
    params:
        th = config['contig_detection']
    run:
        dictlt = {}
        with gzip.open(input.f,'rt') as infile:
            for line in infile:
                lst = line.strip().split()
                if float(lst[1]) > float(params.th):
                    dictlt[lst[0]] = lst[0]

        with gzip.open(input.f2,'rt') as infile:
            with open(output.f,'w') as outfile:
                for i in range(9):
                    infile.readline()
                for line in infile:
                    lst = line.strip().split()
                    try:
                        st = "{}\t{}\n".format(dictlt[lst[0]],lst[1])
                        outfile.write(st)
                    except:
                        st = "{}\t{}\n".format(lst[0],0)
                        outfile.write(st)


rule linear_paste_relative_circu_contig_abundance:
    input:
        f = expand("circular_non_redundant_contig/contig_abundance/per_sample/{sample}_existing_contigs_ab_per_sample.txt",sample=config["samples"])
    output:
        f = "circular_non_redundant_contig/contig_abundance/all_samples_contig_abundance.txt"
    threads: config['threads']
    run:
        if os.path.exists(output.f):
            pass
        else:
            dict = defaultdict(str)
            header = "gene_ID"

            for file in input.f:
                sample_ID = os.path.basename(file)[:-len("_existing_contigs_ab_per_sample.txt")]

            header = "contig_ID"

            for file in input.f:
                sample_ID = os.path.basename(file)[:-len("_existing_contigs_ab_per_sample.txt")]

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

# classify
rule circular_run_mob_typer:
    input:
        f1="circular_non_redundant_contig/circular_non_redundant_contigs"
    output:
        f2="circular_non_redundant_contig/plasmid_classify/plasmid_classify.txt"
    threads: config["threads_mob-typer"]
    conda:
        "%s/mobtyper.yaml" % CONDAENV
    log:
        err="log/circular_run_classify_plasmid/run_mob_typer.err",
        out="log/circular_run_classify_plasmid/run_mob_typer.out"
    shell:
        "mob_typer --multi -n {threads} --infile {input.f1}_rep_seq.fasta --out_file {output.f2} 2>{log.err} >{log.out}"

# rule circular_split_contigs:
#     input:
#         f1="circular_non_redundant_contig/circular_non_redundant_contigs"
#     output:
#         f2=temp(directory("circular_plasmid_contig_split"))
#     threads: config["threads"]
#     run:
#         if not os.path.exists(output.f2):
#             os.mkdir(output.f2)
#         with open(input.f1+"_rep_seq.fasta",mode='r') as infile:
#             dir = defaultdict(str)
#             for line in infile:
#                 if line.startswith(">"):
#                     k = line
#                     dir[k] = ""
#                 else:
#                     dir[k] += line
#
#             for k, v in dir.items():
#                 res = k[1:].split()[0]
#                 a = os.path.join(output.f2, res)+".fna"
#                 with open(a,mode="w") as outfile:
#                     outfile.write(k)
#                     outfile.write(v)
#
#
# rule circular_run_mob_typer:
#     input:
#         f1="circular_plasmid_contig_split"
#     output:
#         f2=temp(directory("circular_non_redundant_contig/all_result_mob_typer"))
#     threads: config["threads_mob-typer"]
#     conda:
#         "%s/mobtyper.yaml" % CONDAENV
#     log:
#         err="log/circular_run_classify_plasmid/run_mob_typer.err",
#         out="log/circular_run_classify_plasmid/run_mob_typer.out"
#     script:
#         "../scripts/run_mob_typer.py"
#
#
#
# rule circular_clean_mob_typer_file:
#     input:
#         f1="circular_non_redundant_contig/all_result_mob_typer"
#     output:
#         f2="circular_non_redundant_contig/plasmid_classify/plasmid_classify.txt"
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



