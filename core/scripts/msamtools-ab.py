#!uer/bin/env python3
# -*- coding:utf-8 -*-

id = os.path.basename(snakemake.input.f)[:-len("_plasmid_sort_filter.bam")]
os.system("msamtools profile --multi=all --unit=" + snakemake.params.ff +
          "  --label=" + id + "  -o " + snakemake.output.f + " " + snakemake.input.f)
