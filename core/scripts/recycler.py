#!uer/bin/env python3
# -*- coding:utf-8 -*-
import os

os.system("samtools index " + snakemake.input[1])
if os.path.exists(snakemake.input[0]):
    os.system("python "+snakemake.params[1]+" -g "+snakemake.input[0]+ " -k "+str(snakemake.params[0])+" -b "+
              snakemake.input[1]+" -o "+snakemake.output[0]+" 2>"snakemake.log.err+" >"+snakemake.log.out)

    os.system("rm -rf "+snakemake.input[1]+".bai")
else:
    print(f"fastg file is not exist, check assembler! ")
