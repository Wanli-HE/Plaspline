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

from lib import common as c


@click.command(
    "working",
    context_settings=dict(ignore_unknown_options=True),
    short_help="run atlas main workflow"
)
@click.argument(
    "workflow",
    default="all",
    type=click.Choice(["qc","assembly","circularized","isolation",
                       "contig_circular","contig_linear",
                       "functional_linear","functional_circular"
                       "gene_circular","gene_linear","all"]),
)
@click.option("-w",
    "--working-dir",
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    help="location to run atlas.",
    default="."
)
@click.option(
    "-j",
    "--jobs",
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help="use at most this many jobs in parallel (see cluster submission for mor details).",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.",
)


@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)


def working(workflow, working_dir, jobs,  dryrun,snakemake_args):
    """
    """

    # c.verify_file()

    cmd = (
        "snakemake --snakefile {snakefile} --directory {working_dir}"
        " --jobs {jobs} --rerun-incomplete"
        " --nolock"
        " --use-conda --conda-prefix {conda_env}"
        " {dryrun}"
        " {target_rule}"
        " {args} "
    ).format(
        snakefile= os.path.join(BASE_PATH,"core","SNAKEMAKE"),
        working_dir=working_dir,
        jobs=jobs,
        dryrun="--dryrun" if dryrun else "",
        args="".join(snakemake_args),
        target_rule=workflow if workflow!="None" else "",
        conda_env= os.path.join(BASE_PATH,'conda_envs')
    )

    logging.info("Executing: %s" % cmd)

    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)

if __name__ == "__main__":
    working()