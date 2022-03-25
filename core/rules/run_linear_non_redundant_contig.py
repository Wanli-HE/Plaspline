#######################################################################################################
# !/usr/bin/env python3
# -*- coding:utf-8 -*-
# pipeline for plasmid and mobile
# 2021-1-21
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


rule non_redundant_nucler_linear_plasmid:
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
        "bwa index {input.f}_rep_seq.fasta 2>{log.err} >{log.out};"
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
        "bwa mem -M {input.f1}_rep_seq.fasta {input.f2} {input.f3} -t {threads} 2>{log.err}|" \
        "samtools view -Sb - >{output.f} "


rule linear_get_sort_file_plasmid:
    input:
        f="linear_non_redundant_contig/{sample}_plasmids.bam"
    output:
        f=temp("linear_non_redundant_contig/{sample}_plasmids_sort.bam")
    threads: config['threads']
    conda:
        "%s/non_redundant.yaml" % CONDAENV
    log:
        err = "log/linear_get_sort_file_plasmid/{sample}.err",
        out = "log/linear_get_sort_file_plasmid/{sample}.out"
    shell:
        "samtools sort -n -@ {threads} {input.f} -o {output.f}  2>{log.err} >{log.out}"


rule linear_reads_mapping:
    input:
        f="linear_non_redundant_contig/{sample}_plasmids_sort.bam"
    output:
        f=temp("linear_non_redundant_contig/{sample}_plasmid_sort_filter.bam")
    threads: config['threads']
    conda:
        "%s/msamtools.yaml" % CONDAENV
    params:
        l = config["linear_msamtools_l"],
        p = config["linear_msamtools_p"],
        z = config["linear_msamtools_z"]
    shell:
        "msamtools filter -b -l {params.l} -p {params.p} -z {params.z} --besthit {input.f} >{output.f}"


rule linear_gene_abundance_msa:
    input:
        f = "linear_non_redundant_contig/{sample}_plasmid_sort_filter.bam"
    output:
        f = temp("linear_non_redundant_contig/{sample}_profile_rb.txt")
    params:
        ff = config["linear_msamtools_contig_abundance_unit"]
    conda:
        "%s/msamtools.yaml" % CONDAENV
    shell:
        "msamtools profile --multi=all --unit={params.ff}  " \
        "--label={input.f} -o {output.f} {input.f}"


rule linear_coverage:
    input:
        f = "linear_non_redundant_contig/{sample}_plasmid_sort_filter.bam"
    output:
        f = temp("linear_non_redundant_contig/{sample}_coverage_info.txt.gz")
    conda:
        "%s/msamtools.yaml" % CONDAENV
    shell:
        "msamtools coverage -z --summary -o {output.f} {input.f}"


rule existing_linear_contigs:
    input:
        f = "linear_non_redundant_contig/{sample}_coverage_info.txt.gz",
        f2 = "linear_non_redundant_contig/{sample}_profile_rb.txt"
    output:
        f = temp("linear_non_redundant_contig/contig_abundance/per_sample/{sample}_existing_contigs_ab_per_sample.txt")
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


rule linear_paste_relative_contig_abundance:
    input:
        f = expand("linear_non_redundant_contig/contig_abundance/per_sample/{sample}_existing_contigs_ab_per_sample.txt",sample=config["samples"])
    output:
        f = "linear_non_redundant_contig/contig_abundance/all_samples_contig_abundance.txt"
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

