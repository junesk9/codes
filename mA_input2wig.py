#!/usr/bin/env python3
"""
Convert Input to Wig

A 5mC *.input file as a reference,
make three *wig files by the seq context (CG, CHG, CHH)

takes ~ 10 min for a file

Usage:
python3 mA_input2wig.py -ch 1,2  -g genome.fa *.input

17th Nov 2014 ToshiMori

=======================Change log

17.12.21 Add a function to parsing restricted chromosomes, comma seperated
14.12.09 Add function to make file 'all'
         a single file contain 5mC level of all sequence contexts
14.12.01 Available various chromosome numbers
         Add option "-g" harbor a multi-fasta genome file,
         when the input file should be corrected (e.g. strand info)
14.11.17 Now only adaptable to A. thaliana, with chromosome# as 5
"""

import sys, time

#######################################################
############# Arguments

args = {"-g":"","-ch":""}

ctime = time.time()
for en,i in enumerate(sys.argv):
    if i == "-g":
        args["-g"] = sys.argv[en+1]
    if i == "-ch":
        args["-ch"] = sys.argv[en+1]


if args["-ch"] != "":
    ch_mode = True
    ch_given = str(args["-ch"])
    ch_given = [i for i in ch_given.split(",")]
else:
    ch_mode = False

fi = sys.argv[-1]
prefix = fi.split("/")[-1]




#######################################################
########### Methods

def parse2wig(dat,ctx,ch):

    name = prefix+"_"+ctx+".wig"
    out = open(name,'w')
    if ctx == "CG":
        col = "color=205,205,0"
    elif ctx == "CHG":
        col = "color=100,149,237"
    elif ctx == "CHH":
        col = "color=205,155,155"
    elif ctx == "all":
        col = "color=255,255,255"
    print("track","type=wiggle_0", "name=%s" % name, "description=%s" % name, "visibility=full", col, "graphType=bar", "viewLimits=-1.0:1.0",sep="\t",file=out)

    ch = ch_given if ch_mode == True else ch ## To restrict the chromosomes output, added 17.12.21
    ch_dic = []
    for i in ch:
        ch_dic.append({})
    for j in dat:
        j[1] = int(j[1])
        idx = ch.index(j[0])
        ch_dic[idx][j[1]] = j[2]

           
    for n,j in enumerate(ch):
        print("variableStep","chrom=%s" % j, sep="\t",file=out)
        seq = sorted(ch_dic[n].keys())
        for k in seq:
            print(k, ch_dic[n][k], sep='\t', file=out)
    out.close()

####################################################
###### Body


### Dilute chr info & Storing input data
ch = []
inp = []
if args["-g"] != "":
    for i in open(args["-g"]):
        if i.startswith(">"):
            i = i.rstrip()[1:]
            ch.append(i)
    ch = sorted(set(ch))
    for j in open(fi):
        chrs = j.split('\t')[0]
        if chrs in ch:
            j = j.strip().split('\t')
            inp.append(j)
else:
    for k in open(fi):
        k = k.strip().split('\t')
        chrs = k[0]
        inp.append(k)
        if not chrs in ch:
            ch.append(chrs)
    ch = sorted(set(ch))
rtime = time.time()
print("data input done, takes %s secs" % int(rtime-ctime))


### Parse by sequence context
CG, CHG, CHH = [],[],[]
allC = []
for i in inp:
    chs = i[0]
    site = i[1]
    ori = i[2]
    rate = 0
    if not int(i[3])+int(i[4]) == 0:
        rate = round(int(i[3])/(int(i[3])+int(i[4])),7)
    if ori == "-":
        rate = rate * -1
    data = [chs,site,rate]
    ctx = i[5]
    if ctx == "CG":
        CG.append(data)
    elif ctx == "CHG":
        CHG.append(data)
    elif ctx == "CHH":
        CHH.append(data)
    allC.append(data)
        
inp = []
rtime = time.time()
print("Data parsing Done, takes %d secs" % int(rtime-ctime))

parse2wig(CG,"CG",ch)
parse2wig(CHH,"CHH",ch)
parse2wig(CHG,"CHG",ch)
parse2wig(allC,"all",ch)
rtime = time.time()

CG, CGH, CHH, allC  = [],[],[],[] ## Make memory empty
print("Make WIGGLE Done, takes %s secs" % int(rtime-ctime))
