#!uer/bin/env python3
# -*- coding:utf-8 -*-
import os
import sys
import logging
from ruamel.yaml import YAML


BASE_PATH = os.path.dirname(os.path.dirname(sys.path[0]))
sys.path.append(BASE_PATH)

CONDAENV = os.path.join(BASE_PATH,"envs")

configfile: os.path.join(BASE_PATH,"temp","config.yaml")

DB_PATH = os.path.join(BASE_PATH,"db") if config["db_dir"] in [None,"none"] else config["db_dir"]

EGGNOG = os.path.join(DB_PATH,'EggNOGV2')
############################################################################################
rule all:
    input:
        f1 = f"{EGGNOG}/eggnog.db",
        f = f"{EGGNOG}/eggnog_proteins.dmnd",
        f2 = os.path.join(DB_PATH,"finished_rgi"),
        f3 = os.path.join(DB_PATH,"VFDB_setA_pro.fas.dmnd"),
        f4 = os.path.join(DB_PATH,"BacMet2_predicted_database.fasta.dmnd"),
        f5 = os.path.join(DB_PATH,"Pfam-A.hmm"),
      #  f6 = os.path.join(DB_PATH,"PPR_Meta_v_1_0"),
        f7 = os.path.join(DB_PATH,"plasmidverify"),
        # f8 = os.path.join(DB_PATH,"recycler"),
        # f9 = os.path.join(DB_PATH,"scapp"),
        # f10 = os.path.join(DB_PATH,"plasclass"),
        f12 = os.path.join(DB_PATH,"db"),
        f13 = os.path.join(DB_PATH,"finished_qc_env"),
        f14 = os.path.join(DB_PATH,"finished_assembly_env"),
        f15 = os.path.join(DB_PATH,"finished_plasforest_env"),
        f16 = os.path.join(DB_PATH,"finished_mobtyper_env"),
        # f17 = os.path.join(DB_PATH,"finished_non_redundant_env"),
        # f18 = os.path.join(DB_PATH,"rfplasmid"),
        # f19 = os.path.join(DB_PATH,"plsdb_.fasta.dmnd"),

        f20 = os.path.join(DB_PATH,"finished_msamtools_env"),
        f21 = os.path.join(DB_PATH,"plasforest"),
        # f22 = os.path.join(DB_PATH,"finished_plasforest_env"),
        # f23 = os.path.join(DB_PATH,"blast")
        f24 = os.path.join(DB_PATH,"non_redundant"),
        f25 = os.path.join(DB_PATH,"finished_scapp_env"),

        #f26 = os.path.join(DB_PATH,"bindash"),
        f27 = os.path.join(DB_PATH,"finished_mmseqs_env"),
        f28 = os.path.join(DB_PATH,"finished_bedtools_env"),
        #f29 = os.path.join(DB_PATH,"DeepVirFinder"),
        #f30 = os.path.join(DB_PATH,"finished_DeepVirFinder_env"),
        f31 = os.path.join(DB_PATH,"finished_plasmidverify_env"),
        f32 = os.path.join(DB_PATH,"finished_platon_env"),
        # f33 = os.path.join(DB_PATH,""),
        f34 = os.path.join(DB_PATH,"rgi")

    run:
        yaml = YAML(typ='safe')
        yaml.default_flow_style = False

        temp_file = os.path.join(BASE_PATH,"conf","default_config.yaml")
        conf_file = os.path.join(BASE_PATH,"temp","config.yaml")

        with open(temp_file,"r") as f:
            conf = yaml.load(f)
        conf["vfdb_database"] = input.f3
        conf["BacMet2_database"] = input.f4
        conf["Pfam_database"] = input.f5
        # conf["ppr-meta_path"] = os.path.join(input.f6,"PPR_Meta")
        conf["plasmidverify_path"] = os.path.join(input.f7,"bin/viralverify")
        # conf["recycler_path"] = os.path.join(input.f8,"bin","recycle.py")
        # conf["scapp_path"] = os.path.join(DB_PATH,"scapp","scapp","scapp.py")
        # conf["plasclass_path"] = os.path.join(input.f10,"classify_fasta.py")
        conf["platondb"] = os.path.join(DB_PATH,"db")
        # conf["rfplasmid_path"] = os.path.join(input.f18,"rfplasmid.py")
        # conf['plsdb_database'] = input.f19
        conf['msamtools_path'] = os.path.join(input.f20,"msamtools-0.9","msamtools")
        conf['plasforest_path'] = input.f21
        # conf['plasforest'] = os.path.join(input.f21)
        # conf['blast'] = os.path.join(input.f23)
        # conf['cdhit-est_path'] = os.path.join(input.f24,"cd-hit-est")
        #conf['bindash_path'] = os.path.join(input.f26,"bindash")
        #conf['deepvirfinder_path'] = input.f29
        conf['rgi_path'] = input.f34

        with open(conf_file,"w") as f1:
            yaml.dump(conf,f1)


rule download_eggNOG_files:
    output:
        f"{EGGNOG}/eggnog.db",
        f"{EGGNOG}/eggnog_proteins.dmnd"
    threads: 1
    conda:
        f"{CONDAENV}/emapper.yaml"
    shell:
        "download_eggnog_data.py -yf --data_dir {EGGNOG}"


rule download_rgi_database:
    output:
        f = os.path.join(DB_PATH,"data")
    threads: 1
    params:
        d= config["rgi_DB_add"]
    shell:
        "wget -O {output.f} {params} --no-check-certificate"


rule tar_rgi_db:
    input:
        f = os.path.join(DB_PATH,"data")
    threads: 1
    output:
        f = os.path.join(DB_PATH,"card.json")
    shell:
        "tar -xvf {input.f} ./card.json;mv ./card.json {output.f}"


rule system_card:
    input:
        f = os.path.join(DB_PATH,"card.json")
    threads: 1
    output:
        f = directory(os.path.join(DB_PATH,"rgi"))
    conda:
        f"{CONDAENV}/rgi.yaml"
    shell:
        "git clone https://github.com/arpcard/rgi;mv rgi {output}"


rule rgidatabase:
    input:
        f = os.path.join(DB_PATH,"rgi"),
        f2 = os.path.join(DB_PATH,"card.json")  
    threads: 1
    output:
        f = temp(touch(os.path.join(DB_PATH,"finished_rgi"))) 
    conda:
        f"{CONDAENV}/rgi.yaml"
    shell:
        "{input.f}/rgi load --card_json {input.f2}"


rule download_gene_aanotation:
    output:
        f1 = os.path.join(DB_PATH,"VFDB_setA_pro.fas"),
        f2 = os.path.join(DB_PATH,"BacMet2_predicted_database.fasta"),
        f3 = os.path.join(DB_PATH,"Pfam-A.hmm")
    threads: 1
    params:
        db_v = config["vfdb_add"],
        db_B = config["BacMet2_add"],
        db_p = config["pfam_add"]
    shell:
        "wget -O ./VFDB_setA_pro.fas.gz {params.db_v};" \
        " gunzip VFDB_setA_pro.fas.gz;" \
        " mv VFDB_setA_pro.fas {output.f1};" \

        "wget -O ./BacMet2_predicted_database.fasta.gz {params.db_B};" \
        "gunzip BacMet2_predicted_database.fasta.gz;" \
        "mv BacMet2_predicted_database.fasta {output.f2};"

        "wget -O ./Pfam-A.hmm.gz {params.db_p}; gunzip Pfam-A.hmm.gz;mv Pfam-A.hmm {output.f3}"


rule makedb_gene_annotation:
    input:
        f1 = os.path.join(DB_PATH,"VFDB_setA_pro.fas"),
        f2 = os.path.join(DB_PATH,"BacMet2_predicted_database.fasta"),
        f3 = os.path.join(DB_PATH,"Pfam-A.hmm")
    threads: 1
    output:
        f1 =os.path.join(DB_PATH,"VFDB_setA_pro.fas.dmnd"),
        f2 =os.path.join(DB_PATH,"BacMet2_predicted_database.fasta.dmnd")
    conda:
        f"{CONDAENV}/non_redundant.yaml"
    shell:
        "diamond makedb --in {input.f1} -d {output.f1};" \
        "diamond makedb --in {input.f2} -d {output.f2}"


# rule d_plasmidverify:
#     output:
#         f = directory(os.path.join(DB_PATH,"plasmidverify"))
#     conda:
#         f"{CONDAENV}/plasmidverify.yaml"
#     params:
#         p=config["plasmidverify_add"]
#     shell:
#         "wget {params.p};unzip master;rm -rf master;" \
#         "mv plasmidVerify-master {output.f}"


rule d_viralverify:
    output:
        f = directory(os.path.join(DB_PATH,"plasmidverify"))
    # conda:
    #     f"{CONDAENV}/plasmidverify.yaml"
    threads: 1
    params:
        p=config["plasmidverify_add"]
    shell:
        "wget {params.p};unzip master;rm -rf master;" \
        "mv viralVerify-master {output.f}"


# rule deepviral:
#     output:
#         f = directory(os.path.join(DB_PATH,"DeepVirFinder"))
#     # conda:
#     #     f"{CONDAENV}/DeepVirFinder.yaml"
#     threads: 1
#     params:
#         p=config["DeepVirFinder_add"]
#     shell:
#         "git clone {params.p};" \
#         "mv DeepVirFinder {output.f}"

# rule envDeepVirFinder:
#     output:
#         f = temp(touch(os.path.join(DB_PATH,"finished_DeepVirFinders_env")))
#     conda:
#         f"{CONDAENV}/DeepVirFinder.yaml"
#     shell:
#         "echo 'finished_DeepVirFinder_env'"

# rule d_recycler:
#     output:
#         f = directory(os.path.join(DB_PATH,"recycler"))
#     conda:
#         f"{CONDAENV}/circularized-recycler.yaml"
#     params:
#         r = config["recycler_add"]
#     shell:
#         "wget {params.r};unzip v0.7;rm -rf v0.7;" \
#         "mv Recycler-0.7 {output.f};"\
#         #"export PYTHONPATH=$PYTHONPAYH:{output.f}"



# rule d_plasclass:
#     output:
#         f = directory(os.path.join(DB_PATH,"plasclass"))
#     conda:
#         f"{CONDAENV}/linearized-plasclass.yaml"
#     params:
#         pc = config["plasclass_add"]
#     shell:
#         "git clone {params.pc};" \
#         "mv PlasClass {output.f}"


# rule d_rfplasmid:
#     output:
#         f = directory(os.path.join(DB_PATH,"rfplasmid")),
#         f2 = directory(os.path.join(DB_PATH,"checkm"))
#     conda:
#         f"{CONDAENV}/linearized-rfplasmid.yaml"
#     params:
#         s = config["rfplasmid_add"],
#         ck = config["checkmdb_add"]
#     shell:
#         "mkdir {output.f2};" \
#         "wget -O chdb {params.ck};" \
#         "mv chdb {output.f2}/chdb;" \
#         "tar zxvf {output.f2}/chdb -C {output.f2};" \
#         "rm -rf {output.f2}/chdb;" \
#         "checkm data setRoot {output.f2};" \
#         "wget {params.s};" \
#         "tar zxvf 0.3;" \
#         "rm -rf 0.3;" \
#         "mv RFPlasmid-0.3 {output.f}"

rule platon_db:
    output:
        f = directory(os.path.join(DB_PATH,"db"))
    # conda:
    #     f"{CONDAENV}/linearized-platon.yaml"
    threads: 1
    params:
        pla = config["platondb_add"]
    shell:
        "wget -O ./db.tar.gz {params.pla};" \
        "tar -xzf db.tar.gz;rm -rf db.tar.gz;" \
        "mv db {output.f}"


#rule bindash_db:
#    output:
#        f = directory(os.path.join(DB_PATH,"bindash"))
#    threads: 1
#    params:
#        bin = config["bindash_add"]
#    shell:
#        "git clone {params.bin};" \
#        "cd bindash;" \
#        "cmake -DCMAKE_BUILD_TYPE=./bindash .;" \
#        "make;" \
#        "cd ../;" \
#        "mv bindash {output.f}"


rule d_plasforest:
    output:
        f = directory(os.path.join(DB_PATH,"plasforest"))
   # conda:
   #     f"{CONDAENV}/linearized-plasforest.yaml"
    params:
        pc = config["plasforest_add"]
    shell:
        "git clone {params.pc};" \
        "mv PlasForest {output.f}; " \
        "cd {output.f}; " \
        "tar -zxvf plasforest.sav.tar.gz; " \
        "bash database_downloader.sh"


rule env_msamtools:
    output:
        f = temp(touch(os.path.join(DB_PATH,"finished_msamtools_env")))
    threads: 1
    conda:
        f"{CONDAENV}/msamtools.yaml"
    shell:
        "echo 'finished_msamtools_env'"

rule env_mmseqs:
    output:
        temp(touch(os.path.join(DB_PATH,"finished_mmseqs_env")))
    conda:
        f"{CONDAENV}/mmseqs2.yaml"
    threads: 1
    shell:
        "echo 'finished_mmseqs2_env'"

rule env_bedtools:
    output:
        temp(touch(os.path.join(DB_PATH,"finished_bedtools_env")))
    conda:
        f"{CONDAENV}/bedtools.yaml"
    threads: 1
    shell:
        "echo 'finished_bedtools_env'"

rule env_qc:
    output:
        temp(touch(os.path.join(DB_PATH,"finished_qc_env")))
    conda:
        f"{CONDAENV}/quality_control.yaml"
    threads: 1
    shell:
        "echo 'finished_qc_env'"

rule env_assembly:
    output:
        temp(touch(os.path.join(DB_PATH,"finished_assembly_env")))
    conda:
        f"{CONDAENV}/assembly.yaml"
    threads: 1
    shell:
        "echo 'finished_assembly_env'"

rule env_plasforest:
    output:
        temp(touch(os.path.join(DB_PATH,"finished_plasforest_env")))
    conda:
        f"{CONDAENV}/linearized-plasforest.yaml"
    threads: 1
    shell:
        "echo 'finished_plasforest_env'"
        
#rule env_plasflow:
#    output:
#        temp(touch(os.path.join(DB_PATH,"finished_linerized_plasflow_env")))
#    conda:
#        f"{CONDAENV}/linearized-plasflow.yaml"
#    threads: 1
#    shell:
#        "echo 'finished_linerized_plasflow_env'"

rule env_mobtyper:
    output:
        temp(touch(os.path.join(DB_PATH,"finished_mobtyper_env")))
    conda:
        f"{CONDAENV}/mobtyper.yaml"
    threads: 1
    shell:
        "echo 'finished_mobtyper_env'"

rule env_non_redundant:
    output:
        f = temp(touch(os.path.join(DB_PATH,"non_redundant")))
    conda:
        f"{CONDAENV}/non_redundant.yaml"
    threads: 1
    shell:
        "echo 'finished_non_redundant_env'"

rule env_scapp:
    output:
        temp(touch(os.path.join(DB_PATH,"finished_scapp_env")))
    conda:
        f"{CONDAENV}/circularized-scapp.yaml"
    threads: 1
    shell:
        "echo 'finished_scapp_env'"

rule env_plasmidverify:
    output:
        f = temp(touch(os.path.join(DB_PATH,"finished_plasmidverify_env")))
    threads: 1
    conda:
        f"{CONDAENV}/plasmidverify.yaml"
    shell:
        "echo 'finished_plasmidverify_env'"

rule env_platon:
    output:
        f = temp(touch(os.path.join(DB_PATH,"finished_platon_env")))
    threads: 1
    conda:
        f"{CONDAENV}/linearized-platon.yaml"
    shell:
        "echo 'finished_platon_env'"


# rule env_DeepVirFinder:
#     output:
#         f = temp(touch(os.path.join(DB_PATH,"finished_DeepVirFinder_env")))
#     threads: 1
#     conda:
#         f"{CONDAENV}/deepvirfinder.yaml"
#     shell:
#         "echo 'finished_DeepVirFinder_env'"








