#!uer/bin/env python3
# -*- coding:utf-8 -*-

import os,glob
from concurrent.futures import ThreadPoolExecutor

if not os.listdir(snakemake.input.f1):
    os.makedirs(snakemake.output.f2)
else:
    def run(file):
            sampleid = os.path.basename(file)[:-len("_split.fa")] + "_mob_typer_result"
            out = os.path.join(snakemake.output.f2, sampleid)
            os.makedirs(out)
            os.system("mob_typer -i " +file+" -n 1" +" -o " +out+ " 2>"+snakemake.log[0]+" >"+snakemake.log[1])

    tf = ThreadPoolExecutor(snakemake.threads)
    files = glob.glob(snakemake.input.f1 + "/*")
    for file in files:
        tf.submit(run,file)
