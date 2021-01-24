#!uer/bin/env python3
# -*- coding:utf-8 -*-

python src/concatenate.py /path/to/catalogue.fna.gz \
       /path/to/assemblies/sample1/contigs.fasta  /path/to/assemblies/sample2/contigs.fasta  [ ... ]


minimap2 -d catalogue.mmi /path/to/catalogue.fna.gz; # make index

minimap2 -t 8 -N 50 -ax sr catalogue.mmi /path/to/reads/sample1.fw.fq.gz /path/to/reads/sample1.rv.fq.gz
    | samtools view -F 3584 -b --threads 8 > /path/to/bam/sample1.bam


vamb --outdir path/to/outdir --fasta /path/to/catalogue.fna.gz
    --bamfiles /path/to/bam/*.bam -o C --minfasta 200000



























