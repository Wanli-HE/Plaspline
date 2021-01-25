#!uer/bin/env python3
# -*- coding:utf-8 -*-
import sys

read_name = "start"
ref_name = "start"

for i in sys.stdin:
    line = i.strip("\n")
    if line.startswith("@"):
        print(line)
    else:
        tmp = line.split("\t")
        if read_name == "start" and ref_name == "start":
            print(line)
            read_name = tmp[0]
            ref_name = tmp[1]
        if read_name != "*":
            if ref_name != tmp[0] or ref_name != tmp[1]:
                print(line)
        else:
            print(line)
        ref_name = tmp[1]
        read_name = tmp[0]