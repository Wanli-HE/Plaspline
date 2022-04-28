#!uer/bin/env python3
# -*- coding:utf-8 -*-

import logging
import multiprocessing
import os
import sys
import subprocess
import click

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_PATH)

from run_preprocessing import preprocessing
from run_core import working
from run_downloading import downloading


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
#@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """
    """
cli.add_command(preprocessing)   #function name
cli.add_command(working)
cli.add_command(downloading)
if __name__ == "__main__":
    cli()
