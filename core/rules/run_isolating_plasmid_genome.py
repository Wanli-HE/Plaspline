#######################################################################################################
# !/usr/bin/env python3
# -*- coding:utf-8 -*-
# pipeline for plasmid and mobile
# 2020-08-21
# wanli he
###########################################################################################################
#
# ruleorder: adding_header>filtering_lesskb>isolating_plasflow>clean_plasflow_result\
#             >isolating_platon>clean_platon_res>predict_plasforest>clean_plasforest\
#             >cut_predict_plasmid>liner_contig_to_plasmid

usetool = config["assembler"]

rule adding_header:
    input:
        f = "assembly_res/{sample}_assembly_res"
    output:
        f = temp("assembly_res/{sample}_contigs_adding_header.fa")
    run:
        sampleid = os.path.basename(input.f)[:-len("_assmebly_res")]

        if usetool == "megahit":
            assemblied = input.f + "/final.contigs.fa"
        elif usetool == "spades":
            assemblied = input.f + "/contigs.fasta"

        with open(assemblied,"r") as infile:
            with open(output.f,"w") as outfile:
                for line in infile:
                    if line.startswith(">"):
                        st = ">"+sampleid+"_"+line[1:]
                        outfile.write(st)
                    else:
                        outfile.write(line)

rule filtering_lesskb:
    input:
        f1= "assembly_res/{sample}_contigs_adding_header.fa"
    output:
        f1=temp(os.path.join(sys.path[0],"{sample}_contigs_1kb.fasta"))
    run:
        with open(input.f1,"r") as infile:
            with open(output.f1,"w") as outfile:
                dict = defaultdict(str)
                for line in infile:
                    if line.startswith(">"):
                        k =line
                        dict[k] = ""
                    else:
                        dict[k] += line

                for k,v in dict.items():
                    lst=v.strip().split()
                    length = 0
                    for ele in lst:
                        length += len(ele)
                    if length >= 1000:
                        outfile.write(k)
                        outfile.write(dict[k])

# rule isolating_plasflow:
#     input:
#         f1="assembly_res/{sample}_contigs_1kb.fasta"
#     output:
#         f1="linear_plasmid_genome/{sample}_res"
#     threads: config['threads']
#     conda:
#         "%s/linearized-plasflow.yaml" % CONDAENV
#     log:
#         out = "log/isolating-linear-contig/{sample}_linear.out",
#         err = "log/isolating-linear-contig/{sample}_linear.err"
#     params:
#         t = config['plasflow_threshold']
#     shell:
#         "PlasFlow.py  --input {input.f1} " \
#                     "--output {output.f1} " \
#                     "--threshold {params.t} " \
#                     "2>{log.err} >{log.out}"

# rule clean_plasflow_result:
#     input:
#         f="linear_plasmid_genome/{sample}_res"
#     output:
#         f1=temp("linear_plasmid_genome/predict_plasmid/{sample}_predict_plasmid_plasflow.fa")
#     run:
#         f = input.f+"_plasmids.fasta"
#         handle = open(f, "rt")
#         f = SeqIO.parse(handle, "fasta")
#         lst = [i for i in f]
#         SeqIO.write(lst, output.f, "fasta")



rule isolating_platon:
    input:
        f = os.path.join(sys.path[0],"{sample}_contigs_1kb.fasta")
    output:
        f = directory("linear_plasmid_genome/{sample}_probs.out")
    threads:
        config['threads']
    conda:
        "%s/linearized-platon.yaml" % CONDAENV
    params:
        mode = config["platon_mode"],
        db = config["platondb"]
    log:
        out = "log/isolating-linear-contig/{sample}_linear.out",
        err = "log/isolating-linear-contig/{sample}_linear.err"
    shell:
        "platon {input.f} " \
                "--db {params.db} " \
                "--mode {params.mode} " \
                "-o {output.f} " \
                "-t {threads} " \
                "-v 2>{log.err} >{log.out}"

# rule clean_platon_res:
#     input:
#         f = "linear_plasmid_genome/{sample}_probs.out",
#         f2 = "assmebly_res/{sample}_contigs_1kb.fasta"
#     output:
#         f = temp("linear_plasmid_genome/predict_plasmid/{sample}_predict_plasmid_platon.fa")
#     run:
#         sampleid = os.path.basename(input.f2)[:-len(".fasta")]
#         f = input.f+"/"+sampleid+".plasmid.fasta"
#
#         handle = open(f, "rt")
#         f = SeqIO.parse(handle, "fasta")
#         lst = [i for i in f]
#         SeqIO.write(lst, output.f, "fasta")



rule predict_plasforest:
    input:
        f = os.path.join(sys.path[0],"{sample}_contigs_1kb.fasta")
    output:
        f = temp(os.path.join(sys.path[0],"{sample}_contigs_1kb.csv"))
    threads: config["threads"]
    params:
        p = config["plasforest_path"]
    conda:
        "%s/linearized-plasforest.yaml" % CONDAENV
    log:
        out = "log/isolating-linear-plasforest/{sample}_linear.out",
        err = "log/isolating-linear-plasforest/{sample}_linear.err"
    #script:
    #    "../scripts/plasforest.py"
    shell:
        "cd {params.p};python3 PlasForest.py -i {input.f} -r -b -f --threads {threads} -o {output.f};cd -"


rule clean_plasforest:
    input:
        f = os.path.join(sys.path[0],"{sample}_contigs_1kb.csv"),
        f2 = os.path.join(sys.path[0],"{sample}_contigs_1kb.fasta")
    output:
        f = "linear_plasmid_genome/{sample}_predict_plasmid_plasforest.fa"
    run:
        lt = []
        with open(input.f,"r") as infile:
            infile.readline()
            for line in infile:
                lst=line.strip().split(",")
                if lst[-1] == "Plasmid":
                    lt.append(lst[0])

        handle = open(input.f2,"rt")
        f = SeqIO.parse(handle,"fasta")
        lst=[i for i in f if i.id in lt]
        SeqIO.write(lst, output.f, "fasta")


rule cut_predict_plasmid:
    input:
        f1 = "linear_plasmid_genome/{sample}_probs.out",
        f2 = os.path.join(sys.path[0],"{sample}_contigs_1kb.fasta"),
        f3 = "linear_plasmid_genome/{sample}_predict_plasmid_plasforest.fa"
    output:
        f = temp("linear_plasmid_genome/{sample}_predict_plasmid.fa")
    run:
        lst = []
        with open(input.f3,"r") as infile:
            for line in infile:
                if line.startswith(">"):
                    lt = line.strip().split()
                    if lt[0][1:] in lst:
                        pass
                    else:
                        lst.append(lt[0][1:])


        sampleid = os.path.basename(input.f2)[:-len(".fasta")]
        f2 = input.f1+"/"+sampleid+".plasmid.fasta"
        with open(f2,"r") as infile:
            for line in infile:
                if line.startswith(">"):
                    lt = line.strip().split()
                    if lt[0][1:] in lst:
                        pass
                    else:
                        lst.append(lt[0][1:])

        handle = open(input.f2, "rt")
        ff= SeqIO.parse(handle, "fasta")
        lsta = [i for i in ff if i.id in lst]
        SeqIO.write(lsta, output.f, "fasta")


# rule removevir:
#     input:
#         f = "linear_plasmid_genome/{sample}_predict_plasmid_temp.fa"
#     output:
#         f = directory("linear_plasmid_genome/{sample}_remov_viral_genome")
#     threads: config['threads']
#     conda:
#         "%s/deepvirfinder.yaml" % CONDAENV
#     log:
#         out = "log/linear_non_redundant_contig/{sample}_removing_linear_vir.out",
#         err = "log/linear_non_redundant_contig/{sample}_removing_linear_vir.err"
#     shell:
#         "python {config[deepvirfinder_path]}/dvf.py " \
#                 "-i {input.f} " \
#                 "-l 1000 " \
#                 "-m {config[deepvirfinder_path]}/models " \
#                 "-c {threads} " \
#                 "-o {output.f}" \
#                 " 2>{log.err} >{log.out}"
#
#
# rule clean_removevir:
#     input:
#         f = "linear_plasmid_genome/{sample}_remov_viral_genome",
#         f2 = "linear_plasmid_genome/{sample}_predict_plasmid_temp.fa"
#     output:
#         f = "linear_plasmid_genome/{sample}_predict_plasmid.fa"
#     params:
#         s = config["deepvirfinder_score"],
#         p = config['deepvirfinder_pvalue']
#     run:
#         file = glob.glob(input.f+"/*")[0]
#         lt = []
#         with open(file,'r') as infile:
#             infile.readline()
#             for line in infile:
#                 lst = line.strip().split()
#                 if float(lst[2]) >= float(params.s) and float(lst[3]) <= float(params.p):
#                     lt.append(lst[0])
#
#         handle = open(input.f2, "rt")
#         ff= SeqIO.parse(handle, "fasta")
#         lsta = [i for i in ff if i.id not in lt]
#         SeqIO.write(lsta, output.f, "fasta")

# rule linear_genecalling_plasmid_gene_distance:
#     input:
#         f = "linear_plasmid_genome/{sample}_predict_plasmid_2.fa"
#     output:
#         f = temp("linear_plasmid_genome/{sample}_nucl_gene.fa")
#     threads: config['threads']
#     conda:
#         "%s/non_redundant.yaml" % CONDAENV
#     params:
#         p = config['g_prodigal_p']
#     log:
#         out = "log/linear_non_redundant_contig/{sample}_genecalling_genecalling.out",
#         err = "log/linear_non_redundant_contig/{sample}_genecalling_genecalling.err"
#     script:
#         "../scripts/gene_calling.py"
#
#
# rule line_plasmid_gene_distance:
#     input:
#         f = "linear_plasmid_genome/{sample}_nucl_gene.fa",
#         f2 = "linear_plasmid_genome/{sample}_predict_plasmid_2.fa"
#     output:
#         f = "linear_plasmid_genome/{sample}_predict_plasmid.fa"
#     run:
#         dictall = defaultdict(list)
#         with open(input.f,'r') as infile:
#             for line in infile:
#                 if line.startswith('>'):
#                     st = line.strip().split()
#                     tt = st[0][1:].rsplit('_',1)
#                     contig =tt[0]
#                     dictall[contig].append(int(st[2]))
#                     dictall[contig].append(int(st[4]))
#
#         dis = {}
#         for k in dictall.keys():
#             dictall[k].sort(reverse=True)
#             dd = []
#             for i in range(len(dictall[k])-1):
#                 d = dictall[k][i] - dictall[k][i+1]
#                 dd.append(d)
#             dis[k] = np.mean(dd)
#
#         nums = [i for i in dis.values()]
#
#         pt = np.percentile(nums, (25, 50, 75), interpolation='midpoint')
#         downp = float(pt[0])
#         upp = float(pt[-1])
#
#         dt = []
#         for k in dis.keys():
#             if upp >= float(dis[k]) and float(dis[k]) >= downp:
#                 dt.append(k)
#
#         handle = open(input.f2, "rt")
#         ff= SeqIO.parse(handle, "fasta")
#         lsta = [i for i in ff if i.id in dt]
#         SeqIO.write(lsta, output.f, "fasta")

rule liner_contig_to_plasmid:
    input:
        f1 = expand("linear_plasmid_genome/{sample}_predict_plasmid.fa",sample=config["samples"]),
        f2 = expand("assembly_res/{sample}_assembly_res",sample=config["samples"])
    output:
        f = "remindering_report/liner_contig_to_plasmid.txt"
    run:
        dict = {}

        for file in input.f1:
            sampleid = os.path.basename(file)[:-len("_predict_plasmid.fa")]

            if usetool == "megahit":
                assemblied = "assembly_res/" + sampleid + "_assembly_res" + "/final.contigs.fa"
            elif usetool == "spades":
                assemblied = "assembly_res/" + sampleid + "_assembly_res" + "/contigs.fasta"


            (contig_num, contig_bp, plasmid_num, plasmid_bp) =(0,0,0,0)

            with open(file,"r") as infile1:
                with open(assemblied,"r") as infile2:
                    for line in infile1:
                        if line.startswith(">"):
                            plasmid_num += 1
                        else:
                            plasmid_bp += len(line.strip())

                    for line in infile2:
                        if line.startswith(">"):
                            contig_num += 1
                        else:
                            contig_bp += len(line.strip())

                    lose = contig_num - plasmid_num
                    rate = plasmid_bp/contig_bp
                    st = f"\t{contig_num}\t{contig_bp}\t{plasmid_num}\t{plasmid_bp}\t{lose}\t{rate}\n"
                    dict[sampleid] = st

        with open(output.f,"w") as outfile:
            header = "sampelid\tassembly_contig_num\tassembly_contig_bp\t" \
                     "plasmid_num\tplasmid_bp\t" \
                     "losing_contig_number\tbase_using_rate\n"
            outfile.write(header)
            for k in dict.keys():
                line = k +dict[k]
                outfile.write(line)
















