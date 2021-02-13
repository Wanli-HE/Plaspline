0#!uer/bin/env python3
# -*- coding:utf-8 -*-

import os,sys


os.mkdir(snakemake.output.f)
out = snakemake.output.f+"/"+os.path.basename(snakemake.input.f)[:-len("_prodigal_protein_seq.faa")] + "_func_annotation_res"
os.system("emapper.py -i "+snakemake.input.f+" --output "+out+" -m "+
          snakemake.params.m+" --cpu "+str(snakemake.threads)+
          " --data_dir "+snakemake.params.d+" 2>"+snakemake.log.err+" >"+snakemake.log.out)
# os.system("cp " + out + "* " + snakemake.output.f)