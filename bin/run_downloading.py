#!uer/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
import subprocess
import click
import logging
import multiprocessing
from ruamel.yaml import YAML

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_PATH)

# Download
@click.command(
    "downloading",
    context_settings=dict(ignore_unknown_options=True),
    short_help="download reference files (need ~50GB)",
)

@click.option(
    "-j",
    "--jobs",
    default=1,
    type=int,
    show_default=True,
    help="number of simultaneous downloads",
)

@click.option(
    "--vfdb",
<<<<<<< HEAD
    help="Virulence Factors Database downloading address:http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz (default) ",
=======
    help="Virulence Factors Database downloading address",
>>>>>>> 76c42201ef2a51b07c231386be73c559bf586e44
    default="http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz"
)

@click.option(
    "--bacmet2",
<<<<<<< HEAD
    help="BacMet Database downloading address: defualt is predicted_database:http:"
         "//bacmet.biomedicine.gu.se/download/BacMet2_predicted_database.fasta.gz (default)",
=======
    help="BacMet Database downloading address: defualt is predicted_database",
>>>>>>> 76c42201ef2a51b07c231386be73c559bf586e44
    default="http://bacmet.biomedicine.gu.se/download/BacMet2_predicted_database.fasta.gz"
)

@click.option(
    "--pfam",
<<<<<<< HEAD
    help="plasmidverify database downloading address:ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/Pfam-A.hmm.gz (default)",
=======
    help="plasmidverify database downloading address",
>>>>>>> 76c42201ef2a51b07c231386be73c559bf586e44
    default="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/Pfam-A.hmm.gz"
)

@click.option(
    "--carddb",
<<<<<<< HEAD
    help="ARGs database address:https://card.mcmaster.ca/latest/data (default)",
=======
    help="ARGs database address",
>>>>>>> 76c42201ef2a51b07c231386be73c559bf586e44
    default="https://card.mcmaster.ca/latest/data"
)

# @click.option(
#     "--pprmeta",
#     help="ppr-meta web address",
#     default="http://cqb.pku.edu.cn/ZhuLab/PPR_Meta/PPR_Meta_v_1_0.zip"
# )

@click.option(
    "--plasmidverify",
<<<<<<< HEAD
    help="plasmidverify  web address:https://codeload.github.com/ablab/plasmidVerify/zip/master (default)",
=======
    help="plasmidverify  web address",
>>>>>>> 76c42201ef2a51b07c231386be73c559bf586e44
    default="https://codeload.github.com/ablab/plasmidVerify/zip/master"
)

# @click.option(
#     "--recycler",
#     help="recycler web address",
#     default="https://codeload.github.com/Shamir-Lab/Recycler/zip/v0.7"
# )

@click.option(
    "--scapp",
<<<<<<< HEAD
    help="scapp web address:https://codeload.github.com/Shamir-Lab/SCAPP/tar.gz/0.1.1 (default)",
=======
    help="scapp web address (v=0.1.1)",
>>>>>>> 76c42201ef2a51b07c231386be73c559bf586e44
    default="https://codeload.github.com/Shamir-Lab/SCAPP/tar.gz/0.1.1"
)

@click.option(
    "--platondb",
<<<<<<< HEAD
    help="platondb web address:https://zenodo.org/record/3924529/files/db.tar.gz (default)",
=======
    help="platondb web address",
>>>>>>> 76c42201ef2a51b07c231386be73c559bf586e44
    default="https://zenodo.org/record/3924529/files/db.tar.gz"
)

# @click.option(
#     "--plasclass",
#     help="plasclass web address",
#     default="https://github.com/Shamir-Lab/PlasClass.git"
# )
#
# @click.option(
#     "--rfplasmid",
#     help="rfplasmid web address (v=0.3)",
#     default="https://codeload.github.com/aldertzomer/RFPlasmid/tar.gz/0.3"
# )

@click.option(
    "--checkmdb",
<<<<<<< HEAD
    help="checkm database:https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz (default)",
=======
    help="checkm database",
>>>>>>> 76c42201ef2a51b07c231386be73c559bf586e44
    default="https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
)

@click.option(
    "--msamtools",
<<<<<<< HEAD
    help="msamtools:https://github.com/arumugamlab/msamtools/releases/download/0.9/msamtools-0.9.tar.gz (default)",
=======
    help="msamtools",
>>>>>>> 76c42201ef2a51b07c231386be73c559bf586e44
    default="https://github.com/arumugamlab/msamtools/releases/download/0.9/msamtools-0.9.tar.gz"
)

@click.option(
    "--plsdb",
<<<<<<< HEAD
    help="plasmid database:https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/?zip (default)",
=======
    help="plasmidome database",
>>>>>>> 76c42201ef2a51b07c231386be73c559bf586e44
    default="https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/?zip"
)

# @click.option(
#     "--plasforest",
#     help="plasforest",
#     default="https://github.com/leaemiliepradier/PlasForest"
# )
#
@click.option(
    "--cdhit",
<<<<<<< HEAD
    help="cdhit:https://github.com/weizhongli/cdhit/releases/download/V4.6.8/cd-hit-v4.6.8-2017-1208-source.tar.gz (default)",
=======
    help="cdhit",
>>>>>>> 76c42201ef2a51b07c231386be73c559bf586e44
    default="https://github.com/weizhongli/cdhit/releases/download/V4.6.8/cd-hit-v4.6.8-2017-1208-source.tar.gz"
)

@click.option(
    "--bindash",
<<<<<<< HEAD
    help="bindash:https://github.com/zhaoxiaofei/bindash.git (default)",
=======
    help="bindash",
>>>>>>> 76c42201ef2a51b07c231386be73c559bf586e44
    default="https://github.com/zhaoxiaofei/bindash.git"
)


@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)

def downloading(jobs,vfdb,bacmet2,pfam,carddb,plasmidverify,
                scapp,platondb,checkmdb,plsdb,msamtools,cdhit,bindash,
                snakemake_args):
    """
    """
    yaml = YAML(typ='safe')
    yaml.default_flow_style = False
    temp_file = os.path.join(BASE_PATH, "conf", "default_config.yaml")
    conf_file = os.path.join(BASE_PATH, "temp", "config.yaml")

    with open(temp_file, "r") as f:
        conf = yaml.load(f)
    conf["vfdb_add"] = vfdb
    conf["BacMet2_add"] = bacmet2
    conf["pfam_add"] = pfam
    conf["rgi_DB_add"] = carddb
    # conf["pprmeta_add"]=pprmeta
    conf["plasmidverify_add"]=plasmidverify
    # conf["recycler_add"] = recycler
    conf["scapp_add"] = scapp
    conf["platondb_add"] = platondb
    # conf["plasclass_add"] = plasclass
    # conf["rfplasmid_add"] = rfplasmid
    conf["checkmdb_add"] = checkmdb
    conf['plsdb_add'] = plsdb
    conf['msamtools_add'] = msamtools
    # conf['plasforest_add'] = plasforest
    # conf['blast_add'] = blast
    conf['cdhit_add'] = cdhit
    conf['bindash_add'] = bindash
    with open(conf_file, "w") as f1:
        yaml.dump(conf, f1)

    cmd = (
        "snakemake --snakefile {snakefile} "
        "--jobs {jobs} --rerun-incomplete "
        "--nolock  --use-conda  --conda-prefix {conda_prefix} "
        "{add_args} "
        "{args}"
    ).format(
        snakefile=os.path.join(BASE_PATH,"core/rules/download.py"),
        jobs=jobs,
        conda_prefix=os.path.join(BASE_PATH,'conda_envs'),
        add_args="" if snakemake_args and snakemake_args[0].startswith("-") else "--",
        args=" ".join(snakemake_args),
    )
    logging.info("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)
if __name__ == "__main__":
    downloading()