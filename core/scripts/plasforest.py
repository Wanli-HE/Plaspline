#!uer/bin/env python3
# -*- coding:utf-8 -*-

import os

infile  = snakemake.input.f
outfile = snakemake.output.f

os.system("cd "+snakemake.params.p+" ; "+"python3 PlasForest.py -i "+
          infile+" -r -b -f --threads "+str(snakemake.threads)+" -o "+outfile+" ;cd -")