#######################################################################################################
# !/usr/bin/env python3
# -*- coding:utf-8 -*-
# pipeline for plasmid and mobile
# 2020-05-19
# wanli he
###########################################################################################################


if config["skip_qc"] in [False,"False"]:

    rule adapter_moving:
        input:
            f1 = lambda wildcards: config["samples"][wildcards.sample]["R1"],
            f2 = lambda wildcards: config["samples"][wildcards.sample]["R2"]
        output:
            f1 = temp("qc_reads/{sample}_adapter_moving_1.fastq.gz"),
            f2 = temp("qc_reads/{sample}_adapter_moving_2.fastq.gz")
        threads: config['threads']
        conda:
            f"{CONDAENV}/quality_control.yaml"
        log:
            out = "log/quality_control/adapter_moving_log/{sample}_adapter_moving.out",
            err = "log/quality_control/adapter_moving_log/{sample}_adapter_moving.err"
        params:
            ktrim = config["adapter_bbduk_r"],
            k = config['adapter_bbduk_k'],
            mink = config['adapter_bbduk_mink'],
            mem = config['adapter_bbduk_mem']
        shell:
            "bbduk.sh -Xmx{params.mem} " \
                    "in={input.f1} " \
                    "in2={input.f2} " \
                    "out1={output.f1} " \
                    "out2={output.f2} " \
                    "ktrim={params.ktrim} " \
                    "k={params.k} " \
                    "mink={params.mink} " \
                    "ref={config[adapter_file]} " \
                    "threads={threads} " \
                    "2> {log.err} > {log.out}"


    rule phix_moving:
        input:
            f1 = "qc_reads/{sample}_adapter_moving_1.fastq.gz",
            f2 = "qc_reads/{sample}_adapter_moving_2.fastq.gz"
        output:
            f1 = temp("qc_reads/{sample}_adapter_phix_moving_1.fastq.gz"),
            f2 = temp("qc_reads/{sample}_adapter_phix_moving_2.fastq.gz")
        threads: config['threads']
        conda:
            f"{CONDAENV}/quality_control.yaml"
        log:
            out = "log/quality_control/phix_moving_log/{sample}_adapter_phix_moving.out",
            err = "log/quality_control/phix_moving_log/{sample}_adapter_phix_moving.err"
        params:
            k = config['phix_k'],
            hdist = config['phix_hdist'],
            mem = config['phix_mem']
        shell:
            "bbduk.sh -Xmx{params.mem} " \
                    "in1={input.f1} " \
                    "in2={input.f2} " \
                    "out1={output.f1} " \
                    "out2={output.f2} " \
                    "k={params.k} " \
                    "hdist={params.hdist} " \
                    "threads={threads} " \
                    "ref={config[phix_file]} " \
                    "2> {log.err} > {log.out}"


    rule quality_window_trimming:
        input:
            f1 = "qc_reads/{sample}_adapter_phix_moving_1.fastq.gz",
            f2 = "qc_reads/{sample}_adapter_phix_moving_2.fastq.gz"
        output:
            f1 = "qc_reads/{sample}_qc_1.fastq.gz",
            f2 = "qc_reads/{sample}_qc_2.fastq.gz",
            f3 = temp("qc_reads/{sample}_adapter_phix_moving_sickle_other_sss.fastq.gz")
        threads: config['threads']
        conda:
            f"{CONDAENV}/quality_control.yaml"
        log:
            out="log/quality_control/quality_window_trimming_log/{sample}_adapter_phix_moving_sickle.out",
            err="log/quality_control/quality_window_trimming_log/{sample}_adapter_phix_moving_sickle.err"
        params:
            outgzip = config['filter_outgzip'],
            t = config['filter_t'],
            l = config['filter_l'],
            mode = config['filter_m']
        shell:
            "sickle  {params.mode} " \
                    "{params.outgzip} " \
                    "-f {input.f1} -r {input.f2} " \
                    "-t {params.t}  " \
                    "-o {output.f1} " \
                    "-p {output.f2} " \
                    "-s {output.f3} " \
                    "-l {params.l} " \
                    "2>{log.err} > {log.out}"

    rule raw_qc_reads:
        input:
            f1 = lambda wildcards: config["samples"][wildcards.sample]["R1"],
            f2 = lambda wildcards: config["samples"][wildcards.sample]["R2"],
            f3 = "qc_reads/{sample}_qc_1.fastq.gz",
            f4 = "qc_reads/{sample}_qc_2.fastq.gz",
        output:
            f = temp("remindering_report/{sample}_raw_reads_to_qc_reads_temp.txt")
        run:
            sampleID = os.path.basename(input.f3)[:-len("_qc_2.fastq.gz")]

            def count(file):
                (reads_number,bp_number) =(0,0)
                try:
                    with gzip.open(file,"r") as infile:
                       index = -1
                       for line in infile:
                           index +=1
                           if index %4 == 0:
                               reads_number += 1
                           elif index %4 == 1:
                               bp_number += float(len(line.strip()))
                except OSError as e:
                    print("check you raw and qc reads fastq file, whether is gzip file!")
                return reads_number,bp_number

            (f_raw_reads,f_raw_bp)=count(input.f1)
            (r_raw_reads,r_raw_bp)=count(input.f2)
            (f_qc_reads,f_qc_bp)=count(input.f3)
            (r_qc_reads,r_qc_bp)=count(input.f4)


            with open(output.f,"w") as outfile:
                st=f"{sampleID}\t{f_raw_reads}\t{f_raw_bp}" \
                    f"\t{r_raw_reads}\t{r_raw_bp}" \
                    f"\t{f_qc_reads}\t{f_qc_bp}" \
                    f"\t{r_qc_reads}\t{r_qc_bp}" \
                    f"\t{f_raw_reads-f_qc_reads}\t{(f_raw_reads-f_qc_reads)/f_raw_reads*100}\t{f_qc_bp/f_raw_bp}" \
                    f"\t{r_raw_reads-r_qc_reads}\t{(r_raw_reads-r_qc_reads)/r_raw_reads*100}\t{r_qc_bp/r_raw_bp}\n"
                outfile.write(st)

    rule catting_raw_qc_reads_adding_header:
        input:
            f =expand("remindering_report/{sample}_raw_reads_to_qc_reads_temp.txt",sample = config["samples"])
        output:
            f = temp("remindering_report/raw_reads_to_qc_reads_temp.txt")
        shell:
            "cat {input.f} >>{output.f}"

    rule raw_qc_reads_adding_header:
        input:
            f = "remindering_report/raw_reads_to_qc_reads_temp.txt"
        output:
            f = "remindering_report/raw_reads_to_qc_reads.txt"
        run:
            with open(input.f,"r") as infile:
                with open(output.f,"w") as outfile:
                    header = "sample_id\tforward_raw_reads_num\tforward_raw_nucleobase_number" \
                             "\treverse_raw_reads_num\treverse_raw_nucleobase_number" \
                             "\tforward_qc_reads_num\tforward_qc_nucleobase_number" \
                             "\treverse_qc_reads_num\treverse_qc_nucleobase_number" \
                             "\tforward_lose_reads_number\tforward_lose_rate_of_reads_num\tforward_using_rate_of_nucleobase" \
                             "\treverse_lose_reads_number\treverse_lose_rate_of_reads_num\treverse_using_rate_of_nucleobase\n"
                    outfile.write(header)
                    for line in infile:
                        outfile.write(line)


else:

    rule changingfilename:
        input:
            f1 = lambda wildcards: config["samples"][wildcards.sample]["R1"],
            f2 = lambda wildcards: config["samples"][wildcards.sample]["R2"]
        output:
            f1 = "qc_reads/{sample}_qc_1.fastq.gz",
            f2 = "qc_reads/{sample}_qc_2.fastq.gz"
        shell:
            "cp {input.f1} {output.f1}; cp {input.f2} {output.f2}"

    rule counting_nothing:
        input:
            f1 = expand("qc_reads/{sample}_qc_1.fastq.gz",sample=config["samples"]),
            f2 = expand("qc_reads/{sample}_qc_2.fastq.gz",sample=config["samples"])
        output:
            f = "remindering_report/raw_reads_to_qc_reads.txt"
        run:
            with open(output.f,"w") as outfile:
                st = "you are skipping quality control step, nothing to collect statistics!!!!"
                outfile.write(st)