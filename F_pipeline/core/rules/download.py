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
        f8 = os.path.join(DB_PATH,"recycler"),
        f9 = os.path.join(DB_PATH,"scapp")
        f10 = os.path.join(DB_PATH,"plasclass"),
        f12 = os.path.join(DB_PATH,"db"),
        f13 = os.path.join(DB_PATH,"finished_qc_env"),
        f14 = os.path.join(DB_PATH,"finished_assembly_env"),
        f15 = os.path.join(DB_PATH,"finished_linerized_plasflow_env"),
        f16 = os.path.join(DB_PATH,"finished_mobtyper_env"),
        f17 = os.path.join(DB_PATH,"finished_non_redundant_env"),
        f18 = os.path.join(DB_PATH,"rfplasmid"),
        f19 = os.path.join(DB_PATH,"plsdb_.fasta.dmnd")
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
        conf["recycler_path"] = os.path.join(input.f8,"bin","recycle.py")
        conf["scapp_path"] = os.path.join(DB_PATH,"scapp","scapp","scapp.py")
        conf["plasclass_path"] = os.path.join(input.f10,"classify_fasta.py")
        conf["platondb"] = os.path.join(DB_PATH,"db")
        conf["rfplasmid_path"] = os.path.join(input.f18,"rfplasmid.py")
        conf['plsdb_database'] = input.f19
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
        "wget -O ./VFDB_setA_pro.fas.gz {params.db_v};" \
        " gunzip VFDB_setA_pro.fas.gz;" \
        " mv VFDB_setA_pro.fas {output.f1};" \

        "wget -O ./BacMet2_predicted_database.fasta.gz {params.db_B};" \
        "gunzip BacMet2_predicted_database.fasta.gz;" \
        "mv BacMet2_predicted_database.fasta {output.f2};"

        "wget -O ./Pfam-A.hmm.gz {params.db_p}; gunzip Pfam-A.hmm.gz; mv Pfam-A.hmm {output.f3};"


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

# rule d_pprmeta:
#     output:
#         f = directory(os.path.join(DB_PATH,"PPR_Meta_v_1_0"))
#     conda:
#         f"{CONDAENV}/linearized-ppr-meta.yaml"
#     params:
#         pm=config["pprmeta_add"]
#     shell:
#         "wget -O ./PPR_Meta_v_1_0.zip {params.pm};" \
#         "unzip PPR_Meta_v_1_0.zip;rm -rf PPR_Meta_v_1_0.zip;" \
#         "mv PPR_Meta_v_1_0 {output.f}"

rule d_recycler:
    output:
        f = directory(os.path.join(DB_PATH,"recycler"))
    conda:
        f"{CONDAENV}/circularized-recycler.yaml"
    params:
        r = config["recycler_add"]
    shell:
        "wget {params.r};unzip v0.7;rm -rf v0.7;" \
        "mv Recycler-0.7 {output.f};"\
        #"export PYTHONPATH=$PYTHONPAYH:{output.f}"

rule d_scapp:
    output:
        f = directory(os.path.join(DB_PATH,"scapp"))
    conda:
        f"{CONDAENV}/circularized-scapp.yaml"
    params:
        s = config["scapp_add"]
    shell:
        "wget {params.s};tar zxvf 0.1.1;rm -rf 0.1.1;" \
        "mv SCAPP-0.1.1 {output.f}"

rule d_plasclass:
    output:
        f = directory(os.path.join(DB_PATH,"plasclass"))
    conda:
        f"{CONDAENV}/linearized-plasclass.yaml"
    params:
        pc = config["plasclass_add"]
    shell:
        "git clone {params.pc};" \
        "mv PlasClass {output.f}"


rule d_rfplasmid:
    output:
        f = directory(os.path.join(DB_PATH,"rfplasmid")),
        f2 = directory(os.path.join(DB_PATH,"checkm"))
    conda:
        f"{CONDAENV}/linearized-rfplasmid.yaml"
    params:
        s = config["rfplasmid_add"],
        ck = config["checkmdb_add"]
    shell:
        "mkdir {output.f2};" \
        "wget -O chdb {params.ck};" \
        "mv chdb {output.f2}/chdb;" \
        "tar zxvf {output.f2}/chdb -C {output.f2};" \
        "rm -rf {output.f2}/chdb;" \
        "checkm data setRoot {output.f2};" \
        "wget {params.s};" \
        "tar zxvf 0.3;" \
        "rm -rf 0.3;" \
        "mv RFPlasmid-0.3 {output.f}"

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

rule plsdb:
    input:
        f = os.path.join(DB_PATH,"db")
    output:
        f = os.path.join(DB_PATH,"plsdb_.fasta.dmnd")
    conda:
        f"{CONDAENV}/non_redundant.yaml"
    params:
        plsdb=config["plsdb_add"]
    shell:
        "wget -O ./plsdb.zip {params.plsdb};" \
        "unzip plsdb.zip;" \
        "blastdbcmd -entry all -db plsdb.fna -out plsdb_.fasta;" \
        "diamond makedb --in plsdb_.fasta -d {output.f};" \
        "rm -rf ./plsdb.* plsdb_changes.tsv README.md plsdb_.fasta"

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
        temp(touch(os.path.join(DB_PATH,"finished_non_redundant_env")))
    conda:
        f"{CONDAENV}/emapper.yaml"
    shell:
        "echo 'finished_emapper_env'"






