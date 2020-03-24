#!/usr/bin/env python3
"""
Parsing "m6_call.DMR.5mc.mp.py" output files
into three files by its sequence context, CG, CHG & CHH
prepared for boxplot() in R

every coloum from single file & Gap filled by 'None'

Oct. 2014 ToshiMori
"""

import glob

fi_ls = glob.glob("*txt")
fi_ls = sorted(fi_ls)
prefix = "_4boxplot.txt"

CG, CHG, CHH = [], [], []
for fi in fi_ls:
    name = fi.split(".")[0]
    column = [[name],[name],[name]]
    print(fi)
    for i in open(fi).readlines()[1:]:
        i = i.rstrip().split("\t")
        column[0].append(i[0])
        column[1].append(i[1])
        column[2].append(i[2])
    CG.append(column[0])
    CHG.append(column[1])
    CHH.append(column[2])

def parsing(ctx, fix):

    name2 = fix+prefix
    out = open(name2,'w')
    ln = []
    for i in ctx:
        while '' in i:
            i.remove('')
        ln.append(len(i))
    max_ln = max(ln)
    ln = len(ctx)

    for i in range(max_ln):
        for j in ctx:
            try: print(j[i], end="\t", file=out)
            except: print("NA", end="\t", file=out)
        print(file=out)
    out.close()

parsing(CG, "CG")
parsing(CHG, "CHG")
parsing(CHH, "CHH")
