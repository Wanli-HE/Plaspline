#######################################################################################################
# !/usr/bin/env python3
# -*- coding:utf-8 -*-
# pipeline for plasmid and mobile
# 2020-08-10
# wanli he
###########################################################################################################

rule get_spliteach_plasmid:
    input:
        f1 = "circular_non_redundant_contig/circular_non_redundant_contigs",
        f2 = config["PLSDB.fasta"]
    output:
        f = temp(directory("circular_non_redundant_contig/all_plasmid_contig_split"))
    run:
        handle = open(input.f1+"_rep_seq.fasta","rt")
        f = SeqIO.parse(handle,"fasta")
        for i in f:
            out = output.f+"/"+i.id+".fasta"
            SeqIO.write(i, out, "fasta")

        handle = open(input.f2, "rt")
        f = SeqIO.parse(handle, "fasta")
        for i in f:
            out = output.f + "/" + i.id + ".fasta"
            SeqIO.write(i, out, "fasta")


rule write_plasmid_each_pathfile:
    input:
        f = temp(directory("circular_non_redundant_contig/all_plasmid_contig_split"))
    output:
        f = temp("circular_non_redundant_contig/plasmid_path.txt")
    run:
        with open(output.f,'w') as outfile:
            files = glob.glob(input.f+"/*")
            for file in files:
                st = file+"\n"
                outfile.write(st)


rule split_plasmid_path_into_groups:
    input:
        f = directory("circular_non_redundant_contig/all_plasmid_contig_split"),
        f2 = "circular_non_redundant_contig/plasmid_path.txt"
    output:
        f = temp(directory("circular_non_redundant_contig/all_plasmid_contig_batches")),
        f2 =temp("circular_non_redundant_contig/batches_comb.txt")
    script:
        "../scripts/split_in_batches.py"


rule pairwise_groups:
    input:
        f = directory("circular_non_redundant_contig/all_plasmid_contig_split"),
        f2 = directory("circular_non_redundant_contig/all_plasmid_contig_batches"),
        f3 = "circular_non_redundant_contig/batches_comb.txt"
    output:
        f = "circular_non_redundant_contig/all_plasmid_host_file.txt"
    conda:
        "%s/mobtyper.yaml" % CONDAENV
    log:
        err="log/circular_run_classify_plasmid/run_mob_typer.err",
        out="log/circular_run_classify_plasmid/run_mob_typer.out"
    script:
        "../scripts/run_mob_typer.py"






