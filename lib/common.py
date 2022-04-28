#!uer/bin/env python3
# -*- coding:utf-8 -*-
# wanli he 0728 2020


import os,sys
import logging
import tempfile
import multiprocessing
from collections import defaultdict

import json
from ruamel.yaml import YAML
from snakemake import utils


BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sample_list_file = "samples.json"
config_file= "config.yaml"

def make_samples_json_file(path,sample_list=sample_list_file):
    """
        creates table sampleID R1 R2 with the absolute paths of fastq files in a given folder
    """
    samples = defaultdict(dict)

    for dir_name, sub_dirs, files in os.walk(os.path.abspath(path)):
        for fname in files:

            if ".fastq" in fname or ".fq" in fname:

                sample_id = fname.split(".fastq")[0].split(".fq")[0]

                sample_id = sample_id.replace("_R1", "").replace("_r1", "").replace("_R2", "").replace("_r2", "")
                # sample_id = sample_id.replace("_1","").replace("_2","")
                sample_id = sample_id.replace("_", "-").replace(" ", "-")

                fq_path = os.path.join(dir_name, fname)

                if "_R2" in fname or "_r2" in fname or "_2" in fname:
                    if 'R2' in samples[sample_id]:
                        logging.error(f"Duplicate sample {sample_id} was found after renaming; "
                                      f"skipping... \n "
                                      f"Samples: \n{samples}")
                    samples[sample_id]['R2'] = fq_path

                elif "_R1" in fname or "_r1" in fname or "_1" in fname:
                    if 'R1' in samples[sample_id]:
                        logging.error(f"Duplicate sample {sample_id} was found after renaming; "
                                      f"skipping... \n "
                                      f"Samples: \n{samples}")
                    samples[sample_id]['R1'] = fq_path

    dit = {"samples":samples}
    f = open(sample_list,'w')
    json.dump(dit,f)
    f.close()
    return sample_list


def make_config_yaml_file(adapter,phix,assembler,skip_qc,config=config_file):
    """
    """

    yaml = YAML(typ='safe')
    yaml.default_flow_style = False

    config_temp = os.path.join(BASE_PATH,"temp","config.yaml")

    with open(config_temp) as config_t:
        conf = yaml.load(config_t)

    # if circularized_tool not in ["metaplasmidspades","recycler", "scapp"]:
    #     logging.critical(f"the param 'circularized_tool'' is {circularized_tool} which is not readable,"
    #                      f"the format should be 'metaplasmidspades','recycler' or 'scapp'!!!!!!")
    #     sys.exit(1)

    if assembler not in ["spades", "megahit"]:
        logging.critical(f"the param 'assembler' is {assembler} which is not readable,"
                         f"the format should be 'spades' or 'megahit' !!!!!!")
        sys.exit(1)

    # if linearized_tool not in ["plasflow", "plasclass", "pprmeta","PFplasmid","platon"]:
    #     logging.critical(f"the param 'linearized_tool'' is {linearized_tool} which is not readable,"
    #                  f"the format should be 'plasflow', 'plasclass', 'cbar' or 'platon' !!!!!!")
    #     sys.exit(1)


    conf["adapter_file"] = adapter
    conf["phix_file"] = phix
    # conf["datatype"] = datatype
    # conf["circularized_tool"] = circularized_tool
    # conf["linearized_tool"] = linearized_tool
    conf["skip_qc"] = skip_qc
    conf["assembler"] = assembler
    if os.path.exists(config):
        logging.warning(f"Config file {config} already exists, I didn't dare to overwrite it. continue...")
    else:
        with open(config, "w") as f:
            yaml.dump(conf, f)
        logging.info(
                     "Configuration file written to %s\n"
                     "You may want to edit it using any text editor."% config
                 )
    return config


def update_config(config_new):

    utils.update_config("config.yaml", config_new)
    return


def verify_file(f1 = "samples.yaml",f2 ="config.yaml", f3 = "Snakefile"):
    BASE_PATH = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    f1_path = 'samples.json'
    f2_path = 'config.yaml'
    f3_path = os.path.join(BASE_PATH,'core','SNAKEMAKE')

    if not os.path.exists(f1_path):
        logging.critical(f"config-file not found: {f1}\n")
                         # "generate one with 'atlas init'")         ######记得改名字
        sys.exit(1)
    if not os.path.exists(f2_path):
        logging.critical(f"config-file not found: {f2}\n")
                         # "generate one with 'atlas init'")
        sys.exit(1)
    if not os.path.exists(f3_path):
        logging.critical(f"config-file not found: {f3}\n")
                         # "generate one with 'atlas init'")
        sys.exit(1)

    # with open(f2_path,"r") as infile:
    #     yaml = YAML(typ='safe')
    #     yaml.default_flow_style = False
    #     conf = yaml.load(infile)
    #     lst=[]
    #     for k in conf:
    #         if conf[k] == None:
    #             lst.append(k)
    #
    #     if len(lst) != 0:
    #         logging.error(f"{[i for i in lst]} is/are still 'None', you need fill these None! in 'temp/config.yaml'")

    return f3_path