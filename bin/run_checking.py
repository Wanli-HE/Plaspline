#!uer/bin/env python3
# -*- coding:utf-8 -*-
import os,glob,sys

a = glob.glob("circular/*_scapp_res")
print(a)
print(len(a))

#path = sys.argv[1]

#files = glob.glob(path)
#print(len(files))
i = 0
for file in a:
    f =  file+"/"+"assembly_graph.confident_cycs.fasta"
    if os.path.exists(f):
        i+=1
    else:
       print(file)
print(i)
