#!/usr/bin/env python3
"""
Common DMC distiller

Find common DMC from TWO DMC files with TWO criteria
1. same chr, sequence
2. the 5mC rates are both positive or negative numbers

The output file is made automatically, with names of two input files.
The corresponding info were adapt from any one of two inputs.
(Look same to other DMC files)

Nov. 13th 2014 ToshiMori

==== Change Log

"""

import sys

##########################################
########## Argument

fi_ls = sys.argv[1:]
name = "comm."+fi_ls[0].split('.')[0]+"."+fi_ls[1].split('.')[0]+".DMC"
out_fi = open(name,'w')
if len(fi_ls) != 2:
    print("""
    #######################################
    Common DMC distiler

    Usage:
    python3 comm.dmc.py DMC1 DMC2


    Nov. 2014 ToshiMori
    #######################################""")
    quit()

###########################################
########### Body
    
fi1 = [i.strip() for i in open(fi_ls[0]).readlines()[1:]]
fi2 = [i.strip().split('\t') for i in open(fi_ls[1]).readlines()[1:]]

out = {}
dic1 = {}
keys1 = []
for i in fi1:
    j = i.split('\t')
    key = j[0]+":"+j[1]
    rate = j[6]
    dic1[key] = float(rate)
    out[key] = i
    keys1.append(key)

dic2 = {}
keys2 = []
for x in fi2:
    key = x[0]+":"+x[1]
    rate = x[6]
    dic2[key] = float(rate)
    keys2.append(key)

keys1 = set(keys1)
keys2 = set(keys2)
comm = keys1.intersection(keys2)
comm = sorted(comm)
print(len(comm))

print("chr","site","s1_mC","s1_NmC","s2_mC","s2_NmC","s2-s1_mC_diff","pValue","context",file=out_fi,sep="\t")
for i in comm:
    if dic1[i] * dic2[i] >  0: ## Alternative use ">" or "<" for select paralell or reversed methylation
        print(out[i],file=out_fi)
out_fi.close()

