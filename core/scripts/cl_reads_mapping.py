#!uer/bin/env python3
# -*- coding:utf-8 -*-

import os

sample_ID = os.path.basename(snakemake.input.f)[:-len("_plasmid_sort.bam")]
out2 = sample_ID + "_names.temp"
out3 = sample_ID + "_count.temp"
os.system("samtools index " + snakemake.input.f)
os.system("samtools idxstats " + snakemake.input.f + "| grep -v '*' | cut -f1 >" + out2)
os.system("samtools idxstats " + snakemake.input.f + "| grep -v '*' | cut -f3 >" + out3)
os.system("paste " + out2 + " " + out3 + " >" + snakemake.output.f)
os.system("rm -rf " + out2 + " " + out3)
temp_file = snakemake.input.f + ".bai"
os.system("rm -rf " + temp_file)