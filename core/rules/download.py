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
        f15 = os.path.join(DB_PATH,"finished_linerized_plasflow_env"),
        f16 = os.path.join(DB_PATH,"finished_mobtyper_env"),
        # f17 = os.path.join(DB_PATH,"finished_non_redundant_env"),
        # f18 = os.path.join(DB_PATH,"rfplasmid"),
        # f19 = os.path.join(DB_PATH,"plsdb_.fasta.dmnd"),
        f20 = os.path.join(DB_PATH,"msamtools"),
        # f21 = os.path.join(DB_PATH,"plasforest"),
        # f22 = os.path.join(DB_PATH,"finished_plasforest_env"),
        # f23 = os.path.join(DB_PATH,"blast")
        f24 = os.path.join(DB_PATH,"cdhit")
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
        conf["plasmidverify_path"] = os.path.join(input.f7,"plasmidverify.py")
        # conf["recycler_path"] = os.path.join(input.f8,"bin","recycle.py")
        # conf["scapp_path"] = os.path.join(DB_PATH,"scapp","scapp","scapp.py")
        # conf["plasclass_path"] = os.path.join(input.f10,"classify_fasta.py")
        conf["platondb"] = os.path.join(DB_PATH,"db")
        # conf["rfplasmid_path"] = os.path.join(input.f18,"rfplasmid.py")
        # conf['plsdb_database'] = input.f19
        conf['msamtools_path'] = os.path.join(input.f20,"msamtools")
        # conf['plasforest_path'] = os.path.join(input.f21,"PlasForest.py")
        # conf['plasforest'] = os.path.join(input.f21)
        # conf['blast'] = os.path.join(input.f23)
        conf['cdhit-est_path'] = os.path.join(input.f24,"cd-hit-est")
        with open(conf_file,"w") as f1:
            yaml.dump(conf,f1)



rule download_eggNOG_files:
    output:
        f"{EGGNOG}/eggnog.db",
        f"{EGGNOG}/eggnog_proteins.dmnd"
    conda:
        f"{CONDAENV}/emapper.yaml"
    shell:
        "download_eggnog_data.py -yf --data_dir {EGGNOG}"


rule download_rgi_database:
    output:
        f = os.path.join(DB_PATH,"data")
    params:
        d= config["rgi_DB_add"]
    shell:
        "wget -O {output.f} {params}"


rule tar_rgi_db:
    input:
        f = os.path.join(DB_PATH,"data")
    output:
        f = os.path.join(DB_PATH,"card.json")
    shell:
        "tar -xvf {input.f} ./card.json;mv ./card.json {output.f}"


rule system_card:
    input:
        f = os.path.join(DB_PATH,"card.json")
    output:
        f = temp(touch(os.path.join(DB_PATH,"finished_rgi")))
    conda:
        f"{CONDAENV}/rgi.yaml"
    shell:
        "rgi load --card_json {input.f}"


rule download_gene_aanotation:
    output:
        f1 = os.path.join(DB_PATH,"VFDB_setA_pro.fas"),
        f2 = os.path.join(DB_PATH,"BacMet2_predicted_database.fasta"),
        f3 = os.path.join(DB_PATH,"Pfam-A.hmm")
    params:
        db_v = config["vfdb_add"],
        db_B = config["BacMet2_add"],
        db_p = config["pfam_add"]
    shell:
        "wget -O ./VFDB_setA_pro.fas.gz {params.db_v} &&" \
        " gunzip VFDB_setA_pro.fas.gz &&" \
        " mv VFDB_setA_pro.fas {output.f1} &&" \

        "wget -O ./BacMet2_predicted_database.fasta.gz {params.db_B}&&" \
        "gunzip BacMet2_predicted_database.fasta.gz&&" \
        "mv BacMet2_predicted_database.fasta {output.f2}&&"

        "wget -O ./Pfam-A.hmm.gz {params.db_p}&& gunzip Pfam-A.hmm.gz&& mv Pfam-A.hmm {output.f3}&&"


rule makedb_gene_annotation:
    input:
        f1 = os.path.join(DB_PATH,"VFDB_setA_pro.fas"),
        f2 = os.path.join(DB_PATH,"BacMet2_predicted_database.fasta"),
        f3 = os.path.join(DB_PATH,"Pfam-A.hmm")
    output:
        f1 =os.path.join(DB_PATH,"VFDB_setA_pro.fas.dmnd"),
        f2 =os.path.join(DB_PATH,"BacMet2_predicted_database.fasta.dmnd")
    conda:
        f"{CONDAENV}/non_redundant.yaml"
    shell:
        "diamond makedb --in {input.f1} -d {output.f1};" \
        "diamond makedb --in {input.f2} -d {output.f2}"


rule d_plasmidverify:
    output:
        f = directory(os.path.join(DB_PATH,"plasmidverify"))
    conda:
        f"{CONDAENV}/plasmidverify.yaml"
    params:
        p=config["plasmidverify_add"]
    shell:
        "wget {params.p};unzip master;rm -rf master;" \
        "mv plasmidVerify-master {output.f}"


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

rule msamtools:
    output:
        f = directory(os.path.join(DB_PATH,"msamtools"))
    conda:
        f"{CONDAENV}/non_redundant.yaml"
    params:
        r = config["msamtools_add"]
    shell:
        "wget -O ./msamtools.tar.gz {params.r} && " \
        "mkdir msamtools && " \
        "tar -zxvf msamtools.tar.gz -C msamtools &&" \
        "cd msamtools/* &&" \
        "./configure;" \
        "make CFLAGS='-std=gnu89 -O2' &&" \
        "cd - && rm -rf masamtools.tar.gz &&" \
        "mv msamtools {output.f}"


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
    conda:
        f"{CONDAENV}/linearized-platon.yaml"
    params:
        pla = config["platondb_add"]
    shell:
        "wget -O ./db.tar.gz {params.pla};" \
        "tar -xzf db.tar.gz;rm -rf db.tar.gz;" \
        "mv db {output.f}"

# rule plsdb:
#     input:
#         f = os.path.join(DB_PATH,"db")
#     output:
#         f = os.path.join(DB_PATH,"plsdb_.fasta.dmnd")
#     conda:
#         f"{CONDAENV}/non_redundant.yaml"
#     params:
#         plsdb=config["plsdb_add"]
#     shell:
#         "wget -O ./plsdb.zip {params.plsdb};" \
#         "unzip plsdb.zip;" \
#         "blastdbcmd -entry all -db plsdb.fna -out plsdb_.fasta;" \
#         "diamond makedb --in plsdb_.fasta -d {output.f};" \
#         "rm -rf ./plsdb.* plsdb_changes.tsv README.md plsdb_.fasta"

rule envqc:
    output:
        temp(touch(os.path.join(DB_PATH,"finished_qc_env")))
    conda:
        f"{CONDAENV}/quality_control.yaml"
    shell:
        "echo 'finished_qc_env'"

rule env_assembly:
    output:
        temp(touch(os.path.join(DB_PATH,"finished_assembly_env")))
    conda:
        f"{CONDAENV}/assembly.yaml"
    shell:
        "echo 'finished_assembly_env'"

rule env_plasflow:
    output:
        temp(touch(os.path.join(DB_PATH,"finished_linerized_plasflow_env")))
    conda:
        f"{CONDAENV}/linearized-plasflow.yaml"
    shell:
        "echo 'finished_linerized_plasflow_env'"

rule env_mobtyper:
    output:
        temp(touch(os.path.join(DB_PATH,"finished_mobtyper_env")))
    conda:
        f"{CONDAENV}/mobtyper.yaml"
    shell:
        "echo 'finished_mobtyper_env'"

rule env_non_redundant:
    output:
        f = directory(os.path.join(DB_PATH,"cdhit"))
    conda:
        f"{CONDAENV}/non_redundant.yaml"
    params:
        cd= config['cdhit_add']
    shell:
        "wget -O ./cdhit.tar.gz {params.cd} && tar -zxvf cdhit.tar.gz && rm -rf cdhit.tar.gz &&" \
        " mv cd-hit-* {output.f} && cd {output.f} && make MAX_SEQ=10000000 &&" \
        " cd cd-hit-auxtools && make MAX_SEQ=10000000"

# rule env_plasforest:
#     input:
#         f =os.path.join(DB_PATH,"blast"),
#         f2 = os.path.join(DB_PATH,"plasforest")
#     output:
#         temp(touch(os.path.join(DB_PATH,"finished_plasforest_env")))
#     conda:
#         f"{CONDAENV}/linearized-plasforest.yaml"
#     shell:
#         "export PATH=$PATH:${input.f}/bin && cd {input.f2} &&" \
#         "tar -zxvf plasforest.sav.tar.gz;bash database_downloader.sh;cd -"
#
# rule download_plasforest:
#     output:
#         f = directory(os.path.join(DB_PATH,"plasforest")),
#         f2 = directory(os.path.join(DB_PATH,"blast"))
#     params:
#         d= config["plasforest_add"],
#         f = config['blast_add']
#     run:
#         a = params.f.rsplit("/",1)[1]
#         b = a.split("+")[0]+"+"
#         os.system("git clone "+params.d+"; mv PlasForest "+output.f+";wget "+params.f+";tar -zxvf "+a+
#                   ";rm -rf "+a+";mv "+b+" "+output.f2)





