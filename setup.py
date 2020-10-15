#!uer/bin/env python3
# -*- coding:utf-8 -*-
from setuptools import setup

setup(
    name='PLASMID PIPELINE', # 应用名
    version="0.0.0", # 版本号
    packages=['F_pipeline'], # 包括在安装包内的 Python 包
)
ruamel.yaml,click,snakemake

# setup(
#     name='metagenome-atlas',
#     version=get_version("atlas/__init__.py"),
#     url='https://github.com/metagenome-atlas/atlas',
#     license='BSD-3',
#     author='Joe Brown, Silas Kieser',
#     author_email='brwnjm@gmail.com, silas.kieser@gmail.com',
#     description='ATLAS - workflows for assembly, annotation, and genomic binning of metagenomic and metatranscriptomic data.',
#     long_description=long_description,
#     long_description_content_type='text/markdown',
#     packages=['atlas'],
#     package_data={'': [
#             "atlas/*",
#                        ]},
#     data_files=[(".", ["README.md", "LICENSE.txt"])],
#     include_package_data=True,
#install_requires= [
# #     ],
#     # install via conda: click, pandas, pyyaml, snakemake
#     entry_points={
#           'console_scripts': [
#               'atlas = atlas.atlas:cli'
#           ]
#     },
#     classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
# )