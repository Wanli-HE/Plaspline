#######################################################################################################
# !/usr/bin/env python3
# -*- coding:utf-8 -*-
# pipeline for plasmid and mobile
# 2020-08-21
# wanli he
###########################################################################################################

usetoolb = config['linearized_tool']


rule adding_header:
    input:
        f = "assmebly_res/{sample}_contigs.fasta"
    output:
        f = temp("assmebly_res/{sample}_contigs_adding_header.fa")
    run:
        sampleid = os.path.basename(input.f)[:-len("_contigs.fasta")]
        with open(input.f,"r") as infile:
            with open(output.f,"w") as outfile:
                for line in infile:
                    if line.startswith(">"):
                        st = ">"+sampleid+"_"+line[1:]
                        outfile.write(st)
                    else:
                        outfile.write(line)


if usetoolb == "plasflow":

    rule filtering_lesskb:
        input:
            f1= "assmebly_res/{sample}_contigs_adding_header.fa"
        output:
            f1=temp("assmebly_res/{sample}_contigs_1kb.fa")
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

    rule isolating_plasflow:
        input:
            f1="assmebly_res/{sample}_contigs_1kb.fa"
        output:
            f1=directory("linear_plasmid_genome/{sample}_res")
        threads: config['threads']
        conda:
            "%s/linearized-plasflow.yaml" % CONDAENV
        log:
            out = "log/isolating-linear-contig/{sample}_linear.out",
            err = "log/isolating-linear-contig/{sample}_linear.err"
        params:
            t = config['plasflow_threshold']
        shell:
            "PlasFlow.py  --input {input.f1} " \
                        "--output {output.f1} " \
                        "--threshold {params.t} " \
                        "2>{log.err} >{log.out}"


    rule clean_plasflow_result:
        input:
            f="linear_plasmid_genome/{sample}_res"
        output:
            f1=temp("linear_plasmid_genome/predict_plasmid/{sample}_predict_plasmid.fa")
        threads: config["threads"]
        shell:
            "cp {input.f}_plasmids.fasta {output.f1}"



elif usetoolb == "plasclass":

    rule isolating_plasclass:
        input:
            f = "assmebly_res/{sample}_contigs_adding_header.fa"
        output:
            f = temp("linear_plasmid_genome/{sample}_probs.out")
        threads: config["threads"]
        conda:
            "%s/linearized-plasclass.yaml" % CONDAENV
        log:
            out = "log/isolating-linear-contig/{sample}_linear.out",
            err = "log/isolating-linear-contig/{sample}_linear.err"
        shell:
            "python {config[plasclass_path]} "
                          "-f {input.f} "
                          "-o {output.f} "
                          "-p {threads} "
                          "2>{log.err} >{log.out}"


    rule isolating_from_prob_out:
        input:
            f1 = "assmebly_res/{sample}_contigs_adding_header.fa",
            f2 = "linear_plasmid_genome/{sample}_probs.out"
        output:
            f1 = temp("linear_plasmid_genome/predict_plasmid/{sample}_predict_plasmid.fa"),
            f2 = temp("linear_plasmid_genome/predict_notplasmid/{sample}_predict_chromosome.fa")
        threads: config["threads"]
        params:
            threshold = config["plasclass_threshold"]
        run:
            t = params.threshold
            dict = defaultdict(str)
            lst_plasmid_contig_id = []

            with open(input.f1,"r") as infile1:
                with open(input.f2,"r") as infile2:
                    with open(output.f1,"w") as outfile1:
                        with open(output.f2,"w") as outfile2:

                            for line in infile1:
                                if line.startswith(">"):
                                    k = line
                                    dict[k] = ""
                                else:
                                    dict[k] += line

                            for line in infile2:
                                lst=line.strip().split()
                                if float(lst[1]) >= float(t):
                                    lst_plasmid_contig_id.append(lst[0])

                            for k in dict.keys():
                                a = k.strip().split()[0][1:]
                                if a in lst_plasmid_contig_id:
                                    st = k+dict[k]
                                    outfile1.write(st)
                                else:
                                    st = k +dict[k]
                                    outfile2.write(st)


elif usetoolb == "platon":
    rule isolating_platon:
        input:
            f = "assmebly_res/{sample}_contigs_adding_header.fa"
        output:
            f = temp(directory("linear_plasmid_genome/{sample}_probs.out"))
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
                    "-c -o {output.f} " \
                    "-t {threads} " \
                    "-v 2>{log.err} >{log.out}"

    rule clean_platon_res:
        input:
            f = "linear_plasmid_genome/{sample}_probs.out"
        output:
            f1 = "linear_plasmid_genome/predict_plasmid/{sample}_predict_plasmid.fa",
            f2 = "linear_plasmid_genome/predict_notplasmid/{sample}_predict_chromosome.fa"
        run:
            sampleid = os.path.basename(input.f)[:-len("_probs.out")]
            f1 = input.f+"/"+sampleid+"_contigs.plasmid.fasta"
            f2 = input.f+"/"+sampleid+"_contigs.chromosome.fasta"
            os.system("cp "+f1+" "+output.f1+";cp "+f2+" "+output.f2)


# elif usetoolb == "rfplasmid":
#     rule copy_folder_and_adding_header:
#         input:
#             f = "assmebly_res"
#         output:
#             f = temp(directory("linear_plasmid_genome/rfplasmid_in"))
#         run:
#             os.makedir(output.f)
#             files = glob.glob(input.f+"/*contigs.fasta")
#             for file in files:
#                 sampleid = os.path.basename(file)[:-len("contigs.fasta")]
#                 out=output.f+"/"+sampleid+"adding_contigs.fasta"
#
#                 with open(file,"r") as infile:
#                     with open(out,"w") as outfile:
#                         for line in infile:
#                             if line.startswith(">"):
#                                 st = ">"+sampleid+line[1:]
#                                 outfile.write(st)
#                             else:
#                                 outfile.write(line)
#
#
#     rule isolating_RFPlasmid:
#         input:
#             f = "linear_plasmid_genome/rfplasmid_in"
#         output:
#             f = "linear_plasmid_genome/rfplasmid_out"
#         threads:
#             config['threads']
#         conda:
#             "%s/linearized-rfplasmid.yaml" % CONDAENV
#         log:
#             out = "log/isolating-linear-contig/linear.out",
#             err = "log/isolating-linear-contig/linear.err"
#         shell:
#             "python3 {config[rfplasmid_path]}   --species {config[rfplamsid_species]} " \
#                                     "--input {input.f} " \
#                                     "--jelly " \
#                                     "--threads {threads} " \
#                                     "--out {output.f} " \
#                                     "2>{log.err} >{log.out}"
#
# elif usetoolb == "ppr-meta":
#     rule isolating_pprmeta:
#             input:
#                 f = "assmebly_res/{sample}_contigs.fasta"
#             output:
#                 f = "linear_plasmid_genome/{sample}_probs.out"
#             threads:
#                 config['threads']
#             conda:
#                 "%s/linearized-ppr-meta.yaml" % CONDAENV
#             params:
#                 t = config["ppr-meta_threshold"]
#             log:
#                 out = "log/isolating-linear-contig/{sample}_linear.out",
#                 err = "log/isolating-linear-contig/{sample}_linear.err"
#             shell:
#                 "{config[ppr-meta_path]} {input.f} " \
#                                     "{output.f} " \
#                                     "-t {params.t} " \
#                                     "2>{log.err} >{log.out}"
#
#     # rule isolating_pprmeta_result:


# rule adding_hearder_precit:
#     input:
#         f = "linear_plasmid_genome/predict_plasmid/{sample}_predict_plasmid.fa"
#     output:
#         f = "linear_plasmid_genome/predict_plasmid/{sample}_predict_plasmid_with_header.fa"
#     run:
#         sampleid = os.path.basename(input.f)[:-len("_predict_plasmid.fa")]
#         with open(input.f,"r") as infile:
#             with open(output.f,"w") as outfile:
#                 for line in infile:
#                     if line.startswith(">"):
#                         st = ">" + sampleid + "_" + line[1:]
#                         outfile.write(st)
#                     else:
#                         outfile.write(line)

rule liner_contig_to_plasmid:
    input:
        f1 = expand("linear_plasmid_genome/predict_plasmid/{sample}_predict_plasmid.fa",sample=config["samples"]),
        f2 = expand("assmebly_res/{sample}_contigs.fasta",sample=config["samples"])
    output:
        f = "remindering_report/liner_contig_to_plasmid.txt"
    run:
        dict = {}
        f_dir = os.path.dirname(input.f2[0])
        for file in input.f1:
            sampleid = os.path.basename(file)[:-len("_predict_plasmid.fa")]
            file2 = f_dir+"/"+sampleid+"_contigs.fasta"
            (contig_num, contig_bp, plasmid_num, plasmid_bp) =(0,0,0,0)

            with open(file,"r") as infile1:
                with open(file2,"r") as infile2:
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







# rule isolation_plasmidverify:
#     input:
#         f1="linear_plasmid_genome/predict_plasmid/{sample}_predict_plasmid.fa"
#     output:
#         f2=temp(directory("linear_plasmid_genome/{sample}_plasmidverify"))
#     threads: config['threads']
#     conda:
#         "%s/plasmidverify.yaml" % CONDAENV
#     params:
#         hmm=config['Pfam_database']
#     log:
#         out = "log/linear_plasmidverify/{sample}_plasmidverify_liner.out",
#         err = "log/linear_plasmidverify/{sample}_plasmidverify_liner.err"
#     shell:
#         "{config[plasmidverify_path]} -f {input.f1} " \
#                              "--hmm {params.hmm} " \
#                              "-t {threads} " \
#                              "-o {output.f2}" \
#                              " 2>{log.err} >{log.out}"
#
#
# rule isolation_moving_chromosome_from_circular:
#     input:
#         f1="linear_plasmid_genome/predict_plasmid/{sample}_predict_plasmid.fa",
#         f2="linear_plasmid_genome/{sample}_plasmidverify"
#     output:
#         f3="linear_plasmid_genome/plasmid/{sample}_verify_plasmid.fasta",
#         f4="linear_plasmid_genome/notplasmid/{sample}_verify_notplasmid.fasta"
#     threads: config['threads']
#     run:
#         file=glob.glob(input.f2+"/*result_table.csv")[0]
#
#         with open(input.f1,mode='r') as infile1:
#             with open(file,mode='r') as infile2:
#                 with open(output.f3,mode='w') as outfile1:
#                     with open(output.f4, mode='w') as outfile2:
#                         lst1=[]
#                         lst2=[]
#                         infile2.readline()
#                         for line in infile2:
#                             lst=line.strip().split(",")
#                             if lst[1] == 'Plasmid':
#                                 lst1.append(lst[0])
#                             elif lst[1] == 'Chromosome':
#                                 lst2.append(lst[0])
#
#                         dict={}
#                         for line in infile1:
#                             if line.startswith(">"):
#                                 k=line
#                                 dict[k]=""
#                             else:
#                                 dict[k] += line
#
#                         for k, v in dict.items():
#                             contigs_ID=k.strip().split()[0][1:]
#                             if contigs_ID in lst1:
#                                 outfile1.write(k)
#                                 outfile1.write(v)
#                             elif contigs_ID in lst2:
#                                 outfile2.write(k)
#                                 outfile2.write(v)


# rule isolation_adding_linear_contig_header:
#     input:
#         f = "linear_plasmid_genome/predict_plasmid/{sample}_predict_plasmid.fa"
#     output:
#         f = "linear_plasmid_genome/predict_plasmid/{sample}_adding_predict_plasmid.fa",
#     threads: config["threads"]
#     run:
#         sampleid = os.path.basename(input.f)[:-len("_predict_plasmid.fa")]
#         with open(input.f,"r") as infile:
#             with open(output.f,"w") as outfile:
#                 for line in infile:
#                     if line.startswith(">"):
#                         st = ">" + sampleid + "_" + line[1:]
#                         outfile.write(st)
#                     else:
#                         outfile.write(line)
#
#         with open(input.f2, "r") as infile:
#             with open(output.f2, "w") as outfile:
#                 for line in infile:
#                     if line.startswith(">"):
#                         st = ">" + sampleid + "_" + line[1:]
#                         outfile.write(st)
#                     else:
#                         outfile.write(line)


#
# rule liner_plasmid_to_verify:
#     input:
#         f1 = expand("linear_plasmid_genome/plasmid/{sample}_verify_plasmid.fasta",sample=config["samples"]),
#         f2 =  expand("linear_plasmid_genome/{sample}_predict_plasmid.fa",sample=config["samples"])
#     output:
#         f = "remindering_report/liner_plasmid_to_verify.txt"
#     run:
#         dict ={}
#         f_dir = os.path.dirname(input.f2[0])
#         for file in input.f1:
#             sampleid = os.path.basename(file)[:-len("_verify_plasmid_linear.fasta")]
#             file2 = f_dir+"/"+sampleid+"_predict_plasmid.fa"
#             (o_lin_num, o_lin_bp, v_lin_num, v_lin_bp) =(0,0,0,0)
#
#             with open(file,"r") as infile1:
#                 with open(file2,"r") as infile2:
#
#                     for line in infile1:
#                         if line.startswith(">"):
#                             v_lin_num += 1
#                         else:
#                             v_lin_bp += len(line.strip())
#
#                     for line in infile2:
#                         if line.startswith(">"):
#                             o_lin_num += 1
#                         else:
#                             o_lin_bp += len(line.strip())
#
#                     lose = o_lin_num - v_lin_num
#                     rate = v_lin_bp/o_lin_bp
#                     st = f"\t{o_lin_num}\t{o_lin_bp}\t{v_lin_num}\t{v_lin_num}\t{lose}\t{rate}\n"
#                     dict[sampleid] = st
#
#         with open(output.f,"w") as outfile:
#             header = "sampelid\tlinear_contig_num(before verifying)\tlinear_contig_bp(before verifying)\t" \
#                      "linear_contig_num(after verifying)\tlinear_contig_bp(after verifying)\t" \
#                      "losing_contig_number\tbase_using_rate(before and after verifying)\n"
#             outfile.write(header)
#             for k in dict.keys():
#                 line = k+dict[k]
#                 outfile.write(line)















