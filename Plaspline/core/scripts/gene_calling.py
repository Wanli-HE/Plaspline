#!uer/bin/env python3
# -*- coding:utf-8 -*-

import os

if os.path.exists(snakemake.input.f):
    os.system("prodigal -p "+snakemake.params.p+" -i "+snakemake.input.f+
              " -d "+snakemake.output.f+" 2>"+snakemake.log.err+" >"+snakemake.log.out)
else:
    pass
