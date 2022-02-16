#!/usr/bin/env python3


import sys
import os.path

suffixes = "conf.dat,gmm.dat,srinfo.dat,counts.dat,gts.dat,expected.dat,paircounts.dat".split(",")


def handleGMM(in_prefix,out_prefix,suffix):
    gmmvf=open(in_prefix+suffix)
    gmmvo=open(out_prefix+suffix,"w")

    head = gmmvf.readline()
    gmmvo.write(head)
    seen_cnvs=set()
    num_line = 0
    for line in gmmvf:
        fields = line.split()
        cnv=fields[0]
        if cnv in seen_cnvs:
            continue
        else:
            seen_cnvs.add(cnv)
        gmmvo.write(line)
        num_line=num_line+1
    gmmvo.close()
    gmmvf.close()

in_prefix=sys.argv[1]
out_prefix=sys.argv[2]

for suffix in suffixes:
    handleGMM(in_prefix,out_prefix,suffix)
