#!uer/bin/env python3
# -*- coding:utf-8 -*-
import os
import sys
import click

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_PATH)


from lib import common as c

@click.command(
    "preprocessing",
    short_help="prepare configuration file and sample table for atlas run",
)

@click.option(
    "--path",
    help = "reads path"
)

@click.option(
    "--adapter",
    help="adapter path"
)

@click.option(
    "--phix",
    help="phix path"
)

@click.option(
    "--datatype",
    type=click.Choice(["metagenome","plasmidome"])
)


@click.option(
    "--circularized_tool",
    type=click.Choice(["metaplasmidspades", "recycler", "scapp"])
)

@click.option(
    "--linearized_tool",
    type=click.Choice(["plasflow", "plasclass", "pprmeta","PFplasmid","platon"])
)

@click.option(
    "--assembler",
    type=click.Choice(["spades", "megahit"])
)

@click.option(
    "--skip_qc",
    default=False,
    help="whether skiping quality control step"
)
# @click.option(
#     "--read_length",
#     default=150
# )
# @click.option(
#     "--config",
#     default="./config.yaml",
#     help="the output path of generated 'config.yaml'"
# )
#
# @click.option(
#     "--sample_list",
#     default="./samples.json",
#     help="the output path of generated 'samples.json'"
# )

def preprocessing(path,datatype,circularized_tool,linearized_tool,adapter,phix,assembler,skip_qc):
    """
    """
    c.make_config_yaml_file(datatype,circularized_tool,linearized_tool,adapter,phix,assembler,skip_qc)
    c.make_samples_json_file(path)

if __name__ == "__main__":
    preprocessing()