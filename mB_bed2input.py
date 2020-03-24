#!/usr/bin/env python3
"""
Convert bed to a 5mC input format file
with revision refer to genome

USAGEe 
python3 mB.bed2input.py bed bed bed ...  genome.fa

Nov.17th 2014 ToshiMori

=================Change log
Feb 16 2015	Also adaptable for Bison output (Three files by sequence context)
Nov 25 2014	Add [try~except] to cope with error in context indexing
		Add filter only accept CG CHG CHH context
Nov 19 2014	minor changes to output file name (including cutoff info)
"""

import glob, sys, time


## Arguments
cutoff = 1 ## Least coverage to identify cytosines to output
ctime = time.time()
fi_ls = sys.argv[1:-1]
gn = sys.argv[-1]
gn = [i for i in open(gn)]

## Process genome.fa 
ch = []
seq = []
k = -1
for i in gn:
    if i.startswith(">"):
        name = i.split(" ")[0].split("\t")[0].split("|")[0]
        name = name[1:].strip()
        ch.append(name)
        seq.append([]) ## In order to convert multi-line sequence to single line
        k += 1
    else:
        i = i.strip()
        seq[k].append(i)

k = 0
for i in seq:
    i = "".join(i)
    seq[k] = i
    k += 1

gn_dic = {}
for n,i in enumerate(ch):
    gn_dic[i] = seq[n]

rtime = time.time()


## Process Bedgraph files 
for fi in fi_ls:
    name = ".".join(fi.split("/")[-1].split(".")[:-1])+".5mC.more"+str(cutoff)+".input"
    out = open(name,'w')
    cnt = "on"
    for i in open(fi).readlines()[1:]:
        i = i.strip().split('\t')
        ch = i[0]
        site = int(i[1])
        mC = int(i[4])
        NmC = int(i[5])
        ## Parse context (ctx) info, whether ctx info. appeared at the filenames(Bison) or including to the bedfile(bwa)
        try:
            ctx = i[6]
        except:
            if cnt == "on":
                print("Warning! no ctx info in coloum6, try to adapt the info. from the filename (Bison-output)")
                cnt = "off"
            ctx = fi.split("/")[-1].split(".")[-2].split("_")[-1]
            #end of if
            if ctx == "CpG":
                ctx = "CG"
            #end of if
        #end of try
        if mC+NmC >= cutoff and ctx in ["CG","CHG","CHH"]:
            ref = gn_dic[ch][site-2:site+3].upper()
            try:
                if ref[2] == "C":
                    ori = "+"
                    ref = ref[2:]
                elif ref[2] == "G":
                    ori = "-"
                    tmp = ref[:3]
                    tmp = tmp[::-1]
                    ref2 = ""
                    for i in tmp:
                        if i == "C":
                            ref2 += "G"
                        elif i == "G":
                            ref2 += "C"
                        elif i == "A":
                            ref2 += "T"
                        elif i == "T":
                            ref2 += "A"
                        else:
                            ref2 += i
                    ref = ref2
		else: 
			sys.exit("Index error no C or G at %s:%s, revise the indexing!!" % (ch,site))
            except:
                print("Warning! error in ctx indexing: %s" % ref)
		if ref[2] == "C"
                	ori = "+"
                elif ref[2] == "G"
			ori = "-"
		ref = "NNN"
		else:
			sys.exit("Index error no C or G at %s:%s, revise the indexing!!" % (ch,site))
            print(ch, site, ori, mC, NmC, ctx, ref, sep="\t", file=out)
            #end of try
        #elif not ctx in ["CG","CHG","CHH"]:
            #print("Warning! cannot recognize the context: %s" % ctx)
        #end of if  
    out.close()  
    rtime = time.time()
    print("Convert %s Done, takes %s secs" % (name, int(rtime-ctime)))
