#!uer/bin/env python3
# -*- coding:utf-8 -*-
import os,sys

os.system("bwa index " + snakemake.input.f1)
f2 = snakemake.input.f1[:-len("_contigs.fasta")] + "_reads_pe.bam"
os.system("bwa mem " + snakemake.input.f1 + " " + snakemake.input.forward + " " + snakemake.input.reverse +
          " -t "+str(snakemake.threads)+" | samtools view -@ "+str(snakemake.threads)+" -buS - > " + f2)
f3 = snakemake.input.f1[:-len("_contigs.fasta")] + "_reads_pe_primary.bam"
os.system("samtools view -@ "+str(snakemake.threads)+" -bF 0x0800 " + f2 + " > " + f3)
os.system("rm -rf " + f2)
os.system("samtools sort -@ "+str(snakemake.threads)+" " + f3 + " -O BAM -o " + snakemake.output.out1)
os.system("rm -rf " + f3)
os.system("rm -rf "_snakemake.input.f1+".*")