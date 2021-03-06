#######################################################################################################
# !/usr/bin/env python3
# -*- coding:utf-8 -*-
# pipeline for plasmid and mobile
# 2020-08-10
# wanli he
###########################################################################################################

usetool = config["assembler"]

if usetool == "megahit":
    rule megahit:
        input:
            f1 = "qc_reads/{sample}_qc_1.fastq.gz",
            f2 = "qc_reads/{sample}_qc_2.fastq.gz"
        output:
            f = directory("assembly_res/{sample}_assembly_res")
        threads: config['threads']
        conda:
            f"{CONDAENV}/assembly.yaml"
        params:
            min_count = config["megahit_min_count"],
            k_min = config["megahit_k_min"],
            k_max = config["megahit_k_max"],
            k_step = config["megahit_k_step"],
            merge_level = config["megahit_merge_level"],
            prune_level = config["megahit_prune_level"],
            low_local_ratio = config["megahit_low_local_ratio"],
            min_contig_len = config["prefilter_minimum_contig_length"],
            preset = config['megahit_preset'],
            m = config['megahit_m']
        log:
            out = "log/assembly_res/{sample}.out",
            err = "log/assembly_res/{sample}.err"
        shell:
            "megahit -1 {input.f1} -2 {input.f2} \
                        --num-cpu-threads {threads} \
                        --k-min {params.k_min} \
                        --k-max {params.k_max} \
                        --k-step {params.k_step} \
                        -o {output.f} \
                        --min-contig-len {params.min_contig_len} \
                        --min-count {params.min_count} \
                        --merge-level {params.merge_level} \
                        --prune-level {params.prune_level} \
                        --low-local-ratio {params.low_local_ratio} \
                        -m {params.m}  \
                        {params.preset}  2>{log.err} >{log.out}"


    # rule cleaning_assembler:
    #     input:
    #         f="assmebly_res/{sample}_assmebly_res"
    #     output:
    #         f="assmebly_res/{sample}_contigs.fasta"
    #     run:
    # file=input.f+"/final.contigs.fa"
    # handle = open(file,"rt")
    # f = SeqIO.parse(handle,"fasta")
    # lst=[i for i in f]
    # SeqIO.write(lst, output.f, "fasta")



else:
    rule SPAdes:
        input:
            f1 = "qc_reads/{sample}_qc_1.fastq.gz",
            f2 = "qc_reads/{sample}_qc_2.fastq.gz"
        output:
            f = directory("assembly_res/{sample}_assembly_res")
        threads: config["threads"]
        conda:
            f"{CONDAENV}/assembly.yaml"
        log:
            out = "log/assembly_res/{sample}.out",
            err = "log/assembly_res/{sample}.err"
        params:
            k = config['spades_k'],
            m = config['spades_m'],
            p = config['spades_preset']
        shell:
            "spades.py  -1 {input.f1} " \
                    "-2 {input.f2} " \
                    "--{params.p} " \
                    "-m {params.m} " \
                    "-t {threads} " \
                    "-k {params.k} " \
                    "-o {output.f} " \
                    "2>{log.err} > {log.out}"

    # rule cleaning_assembler:
    #     input:
    #         f="assmebly_res/{sample}_spades"
    #     output:
    #         f = temp("assmebly_res/{sample}_contigs.fasta"),
    #         f2= temp("assmebly_res/{sample}_contigs.fastg")
    #     run:
    #         f1 = input.f+"/contigs.fasta"
    #         f2 = input.f+"/assembly_graph.fastg"
    #
    #         handle1 = open(f1, "rt")
    #         f1 = SeqIO.parse(handle1, "fasta")
    #         lst1 = [i for i in f1]
    #         SeqIO.write(lst1, output.f, "fasta")
    #
    #         with open(f2,"r") as infile:
    #              with open(output.f2,"w") as outfile:
    #                  for line in infile:
    #                      outfile.write(line)

rule qc_to_assembly:
    input:
        f1 = expand("qc_reads/{sample}_qc_1.fastq.gz",sample=config["samples"]),
        f2 = expand("qc_reads/{sample}_qc_2.fastq.gz",sample=config["samples"]),
        f3 = expand("assembly_res/{sample}_assembly_res",sample=config["samples"])
    output:
        f = "remindering_report/qc_to_assembly.txt"
    run:
        dict = {}
        for file in input.f1:
            sampleid = os.path.basename(file)[:-len("_qc_1.fastq.gz")]
            reverse = file[:-len("_qc_1.fastq.gz")]+"_qc_2.fastq.gz"

            if usetool == "megahit":
                assemblied = "assembly_res/"+sampleid+"_assembly_res"+"/final.contigs.fa"
            elif usetool == "spades":
                assemblied = "assembly_res/"+sampleid+"_assembly_res" + "/contigs.fasta"

            (f_reads_bp,r_reads_bp,contig_bp) = (0,0,0)

            with gzip.open(file,"r") as infile1:
                with gzip.open(reverse,"r") as infile2:
                    with open(assemblied, "r") as infile3:
                        (findex,rindex) = (-1,-1)
                        for line in infile1:
                            findex += 1
                            if findex %4 == 1:
                                f_reads_bp += len(line.strip())

                        for line in infile2:
                            rindex += 1
                            if rindex %4 == 1:
                                f_reads_bp += len(line.strip())

                        for line in infile3:
                            if line.startswith(">"):
                                pass
                            else:
                                contig_bp += len(line.strip())

                        reads_bp = f_reads_bp + r_reads_bp
                        usingrate= contig_bp/reads_bp
                        st = f"\t{reads_bp}\t{contig_bp}\t{usingrate}\n"
                        dict[sampleid] = st

        with open(output.f,"w") as outfile:
            header = "sampleid\treads_bp\tcontig_bp\tusing_rate_of_bp\n"
            outfile.write(header)
            for k in dict.keys():
                line = k+dict[k]
                outfile.write(line)

