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
    help = "reads path,raw reads or QC reads"
)

@click.option(
    "--adapter",
    help="adapter path,if skip Quality Control step, not necessary to input",
    default="None"
)

@click.option(
    "--phix",
    help="phix path,if skip Quality Control step, not necessary to input",
    default="None"
)

# @click.option(
#     "--datatype",
#     type=click.Choice(["metagenome","plasmidome"])
# )
#
#
# @click.option(
#     "--circularized_tool",
#     default="all",
#     type=click.Choice(["metaplasmidspades", "scapp","all"])
# )
#
# @click.option(
#     "--linearized_tool",
#     default="all",
#     type=click.Choice(["plasflow", "plasforest", "platon","all"])
# )

@click.option(
    "--assembler",
    type=click.Choice(["spades", "megahit"]),
    default="spades"
)

@click.option(
    "--skip_qc",
    help="whether skiping quality control step",
    default=False
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

def preprocessing(path,adapter,phix,assembler,skip_qc):
    """
    """
    c.make_config_yaml_file(adapter,phix,assembler,skip_qc)
    c.make_samples_json_file(path)

if __name__ == "__main__":
    preprocessing()