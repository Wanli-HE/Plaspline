#!uer/bin/env python3
# -*- coding:utf-8 -*-

import os,sys


os.mkdir(snakemake.output.f)
sample = os.path.basename(snakemake.input.f)[:-len("_prodigal_protein_seq.faa")] + "_res"
os.system("emapper.py -i "+snakemake.input.f+" --output "+sample+" -m "+
          snakemake.params.m+" --cpu "+str(snakemake.threads)+
          " --data_dir "+snakemake.params.d+" 2>"+snakemake.log.err+" >"+snakemake.log.out)
os.system("mv " + sample + "* " + snakemake.output.f)