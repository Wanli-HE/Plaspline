#######################################################################################################
# !/usr/bin/env python3
# -*- coding:utf-8 -*-
# pipeline for plasmid and mobile
# 2019-03-06
# wanli he
###########################################################################################################


rule metaplasmidSPAdes:
    input:
        f1 = "qc_reads/{sample}_qc_1.fastq.gz",
        f2 = "qc_reads/{sample}_qc_2.fastq.gz"
    output:
        f = directory("circular/{sample}_metaplasmidspades")
    threads: config["threads"]
    conda:
        "%s/assembly.yaml" % CONDAENV
    log:
        out="log/metaplasmidSPAdes/{sample}_metaplasmideSPAdes.out",
        err="log/metaplasmidSPAdes/{sample}_metaplasmideSPAdes.err"
    params:
        k = config['c_spades_k'],
        m = config['c_spades_m'],
        preset1 = config['c_sapdes_preset1'],
        preset2 = config['c_sapdes_preset2']
    shell:
        "spades.py -1 {input.f1} " \
                "-2 {input.f2} " \
                "-m {params.m} " \
                "-t {threads} " \
                "-k {params.k} " \
                "--{params.preset1} " \
                "--{params.preset2} " \
                "-o {output.f} " \
                "2>{log.err} > {log.out}"

# rule cleaning_metaplasmidspades:
#     input:
#         f="{sample}_metaplasmidspades"
#     output:
#         f=temp("circular/{sample}_cycs_metaplasmidspades.fasta")
#     threads: config["threads"]
#     shell:
#         "mv {input.f}/contigs.fasta {output.f}"


if config["assembler"] == "megahit":
    rule make_fastg:
        input:
            f="assembly_res/{sample}_assembly_res"
        output:
            f=temp("assembly_res/{sample}_contigs.fastg")
        threads: config['threads']
        conda:
            "%s/assembly.yaml" % CONDAENV
        params:
            len = config["megahit_core_len"]
        log:
            out = "log/megahit_core/{sample}_megahit.out",
            err = "log/megahit_core/{sample}_megahit.err"
        shell:
            "megahit_core contig2fastg {params.len} {input.f}/final.contigs.fa > {output.f}"

    rule scapp:
        input:
            f1 = "assembly_res/{sample}_contigs.fastg",
            f2 = "qc_reads/{sample}_qc_1.fastq.gz",
            f3 = "qc_reads/{sample}_qc_2.fastq.gz"
        output:
            f1= directory("circular/{sample}_scapp_res")
        threads: config["threads"]
        conda:
            "%s/circularized-scapp.yaml" % CONDAENV
        params:
            k = config["scapp_k"]
        log:
            out=os.path.join("log","scapp","{sample}_scapp.out"),
            err=os.path.join("log","scapp","{sample}_scapp.err")
        shell:
            "scapp -g {input.f1} " \
                    "-o {output.f1} " \
                    "-k {params.k} " \
                    "-r1 {input.f2} " \
                    "-r2 {input.f3} " \
                    "-p {threads} " \
                    "2>{log.err} > {log.out}"

elif config["assembler"] == "spades":
    rule filer_fastg:
        input:
            f = "assembly_res/{sample}_assembly_res"
        output:
            f = temp("assembly_res/{sample}_filter.fastg")
        run:
            f = input.f+"/assembly_graph.fastg"
            handle = open(f, "rt")
            f = SeqIO.parse(handle, "fasta")
            lst = [i for i in f if len(i.seq)>1000]
            SeqIO.write(lst, output.f, "fasta")

    rule scapp:
        input:
            f1 = "assembly_res/{sample}_filter.fastg",
            f2 = "qc_reads/{sample}_qc_1.fastq.gz",
            f3 = "qc_reads/{sample}_qc_2.fastq.gz"
        output:
            f1= directory("circular/{sample}_scapp_res")
        threads: config["threads"]
        conda:
            "%s/circularized-scapp.yaml" % CONDAENV
        params:
            k = config["scapp_k"]
        log:
            out=os.path.join("log","scapp","{sample}_scapp.out"),
            err=os.path.join("log","scapp","{sample}_scapp.err")
        shell:
            "scapp -g {input.f1} " \
                    "-o {output.f1} " \
                    "-k {params.k} " \
                    "-r1 {input.f2} " \
                    "-r2 {input.f3} " \
                    "-p {threads} " \
                    "2>{log.err} > {log.out}"


# rule clean_scapp_res:
#     input:
#         f1="circular/{sample}_scapp_res"
#     output:
#         f1=temp("circular/{sample}_cycs_scapp.fasta")
#     threads: config["threads"]
#     run:
#         folder_name=os.path.basename(input.f1)[:-len("_scapp_res")]
#         file_cir=folder_name+"_contigs.confident_cycs.fasta"
#         os.system("cp "+input.f1+"/"+file_cir+" "+output.f1)
#         os.system("rm -rf "+input.f1)


rule cat_metaplasmidspades_scapp:
    input:
        f1 = "circular/{sample}_metaplasmidspades",
        f2 = "circular/{sample}_scapp_res"
    output:
        f = temp("circular/{sample}_cycs.fasta")

    shell:
        "cat {input.f1}/contigs.fasta {input.f2}/assembly_graph.confident_cycs.fasta > {output.f}"
        #"cat {input.f1}/contigs.fasta {input.f2}/{wildcards.sample}_contigs.confident_cycs.fasta > {output.f}"

    # run:
    #     file1 = input.f1 + "/contigs.fasta"
    #
    #     f = os.path.basename(input.f2)[:-len("_scapp_res")]
    #     file2 = input.f2+"/"+f+"_contigs.confident_cycs.fasta"
    #
    #     os.system("cat "+file1+" "+file2+" >"+output.f)

rule qc_to_circualr_metaplasmidspades:
    input:
        f1 = expand("qc_reads/{sample}_qc_1.fastq.gz",sample=config["samples"]),
        f2 = expand("qc_reads/{sample}_qc_2.fastq.gz",sample=config["samples"]),
        f3 = expand("circular/{sample}_cycs.fasta",sample=config["samples"])
    output:
        f = "remindering_report/qc_to_circualr_metaplasmidspades.txt"
    run:
        dict = {}
        f_dir = os.path.dirname(input.f3[0])
        for file in input.f1:
            sampleid = os.path.basename(file)[:-len("_qc_1.fastq.gz")]
            reverse = file[:-len("_qc_1.fastq.gz")]+"_qc_2.fastq.gz"
            assemblied = f_dir+"/"+sampleid+"_cycs.fasta"
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
                line = k +dict[k]
                outfile.write(line)



rule plasmidverify:
    input:
        f1="circular/{sample}_cycs.fasta"
    output:
        f2=temp(directory("circular/{sample}_plasmidverify"))
    threads: config['threads']
    conda:
        "%s/plasmidverify.yaml" % CONDAENV
    params:
        hmm=config['Pfam_database'],
        threshold = config['verify_threshold']
    log:
        out = "log/circular_plasmidverify/{sample}_plasmidverify_liner.out",
        err = "log/circular_plasmidverify/{sample}_plasmidverify_liner.err"
    shell:
        "{config[plasmidverify_path]} -f {input.f1} " \
                             "--hmm {params.hmm} " \
                             "-t {threads} --thr {params.threshold} " \
                             "-o {output.f2}" \
                             " 2>{log.err} >{log.out}"


rule moving_chromosome_from_circular:
    input:
        f1="circular/{sample}_cycs.fasta",
        f2="circular/{sample}_plasmidverify"
    output:
        f3=temp("circular/plasmid/{sample}_verify_circular_plasmid.fasta"),
        f4=temp("circular/chromosome/{sample}_verify_not_circular_plasmid.fasta")
    threads: config['threads']
    run:
        file=glob.glob(input.f2+"/*result_table.csv")[0]

        with open(input.f1,mode='r') as infile1:
            with open(file,mode='r') as infile2:
                with open(output.f3,mode='w') as outfile1:
                    with open(output.f4, mode='w') as outfile2:
                        lst1=[]
                        lst2=[]
                        infile2.readline()
                        for line in infile2:
                            lst=line.strip().split(",")
                            if lst[1] == 'Plasmid':
                                lst1.append(lst[0])
                            elif lst[1] == 'Chromosome':
                                lst2.append(lst[0])

                        dict={}
                        for line in infile1:
                            if line.startswith(">"):
                                k=line
                                dict[k]=""
                            else:
                                dict[k] += line

                        for k, v in dict.items():
                            contigs_ID=k.strip().split()[0][1:]
                            if contigs_ID in lst1:
                                outfile1.write(k)
                                outfile1.write(v)
                            elif contigs_ID in lst2:
                                outfile2.write(k)
                                outfile2.write(v)


rule adding_circular_contig_header:
    input:
        f = "circular/plasmid/{sample}_verify_circular_plasmid.fasta",
        f2 = "circular/chromosome/{sample}_verify_not_circular_plasmid.fasta"
    output:
        f = "circular/plasmid/{sample}_verify_plasmid_circular.fasta",
        f2 = "circular/chromosome/{sample}_verify_chromosome_circular.fasta"
    threads: config["threads"]
    run:
        sampleid = os.path.basename(input.f)[:-len("_verify_circular_plasmid.fasta")]
        with open(input.f,"r") as infile:
            with open(output.f,"w") as outfile:
                for line in infile:
                    if line.startswith(">"):
                        st = ">" + sampleid + "_" + line[1:]
                        outfile.write(st)
                    else:
                        outfile.write(line)

        with open(input.f2, "r") as infile:
            with open(output.f2, "w") as outfile:
                for line in infile:
                    if line.startswith(">"):
                        st = ">" + sampleid + "_" + line[1:]
                        outfile.write(st)
                    else:
                        outfile.write(line)

rule circular_to_verify:
    input:
        f1 = expand("circular/plasmid/{sample}_verify_circular_plasmid.fasta",sample=config["samples"]),
        f2 = expand("circular/{sample}_cycs.fasta",sample=config["samples"])
    output:
        f = "remindering_report/circualr_to_verify.txt"
    run:
        dict = {}
        f_dir = os.path.dirname(input.f2[0])
        for file in input.f1:
            sampleid = os.path.basename(file)[:-len("_verify_circular_plasmid.fasta")]
            file2 = f_dir+"/"+sampleid+"_cycs.fasta"
            (o_cir_num, o_cir_bp, v_cir_num, v_cir_bp) =(0,0,0,0)

            with open(file,"r") as infile1:
                with open(file2,"r") as infile2:
                    for line in infile1:
                        if line.startswith(">"):
                            v_cir_num += 1
                        else:
                            v_cir_bp += len(line.strip())

                    for line in infile2:
                        if line.startswith(">"):
                            o_cir_num += 1
                        else:
                            o_cir_bp += len(line.strip())

                    lose = o_cir_num - v_cir_num
                    rate = v_cir_bp/o_cir_bp
                    st = f"\t{o_cir_num}\t{o_cir_bp}\t{v_cir_num}\t{v_cir_num}\t{lose}\t{rate}\n"
                    dict[sampleid] = st

        with open(output.f,"w") as outfile:
            header = "sampelid\tlinear_contig_num(before verifying)\tlinear_contig_bp(before verifying)\t" \
                     "linear_contig_num(after verifying)\tlinear_contig_bp(after verifying)\t" \
                     "losing_contig_number\tbase_using_rate(before and after verifying)\n"
            outfile.write(header)
            for k in dict.keys():
                line = k+dict[k]
                outfile.write(line)

