#!/usr/bin/env python3
"""
small RNA counter with Null-model replicator
Now only for acceptable Arabidopsis genome (5 chromosomes)
Two function:

Dec.13th. 2014 June
"""


import random, sys

##############################################################
######## ARGS


N_rep = 1000

arg = sys.argv[1:]
loc_fi = arg[-1]  ## intact *.sam file

args = {"-dmr":""}
for en,i in enumerate(arg):
    if i == "-dmr":
        args["-dmr"] = arg[en+1]

prefix1 = loc_fi.split("/")[-1].split(".")[0]
prefix2 = "random"
if args["-dmr"] != "":
    prefix2 = ".".join(args["-dmr"].split(".")[:-1])
out_fi = prefix1+"."+prefix2+"-NullModelTest_"+str(N_rep)+"rep.txt"

##############################################################
############# Method

def MakeOverlap(group,len_, ctx):

    count = 0
    for idx,i in enumerate(group):
        for j in i:
            for k in at_rng[idx]:
                comm = j.intersection(k)
                if len(comm) > 0:
                    count += 1
    countPer = count /len_ * 1000
    print("Calculate %s done" % ctx)
    print("smRNA overlap to %s DMR: %s" % (ctx, countPer), file=out_fi)
          

##############################################################
####### body

loc = open(loc_fi)
at = {1:30427671, 2:19698289, 3:23459830, 4:18585056, 5:18585056}
at_ch = sorted([i for i in at.keys()])
at_rng = [[],[],[],[],[]]
for i in loc:
    if i.startswith("@"):
        pass
    else:
        i = i.strip().split("\t")
        if not i[2] == "*":
            ch = int(i[2])
            st = int(i[3])
            ed = st+len(i[9])+1
            rng = set(range(st,ed))
            idx = at_ch.index(ch)
            at_rng[idx].append(rng)
print("load SAM file done!")


at = {1:30427671, 2:19698289, 3:23459830, 4:18585056, 5:18585056}

if args["-dmr"] == "":
    out_fi = open(out_fi,'w')
    k= 0
    x,y = 0,100
    result = []
    while k < N_rep:
        ran_ch = random.randrange(len(at))+1
        ran_seq = random.randrange(at[ran_ch])+1
        if ran_seq > 1000:
            ran_rng = set(range(ran_seq-1000,ran_seq))
            idx = at_ch.index(ran_ch)
            count = 0
            for l in at_rng[idx]:
                comm = l.intersection(ran_rng)
                if len(comm) != 0:
                    count += 1
            result.append(count)
            k += 1
            x += 1
        else:
            pass
        if x > 99:
            print("Null Model test, %s reps done!" % y)
            y += 100
            x = 0
    for m in result:
        print(str(m), file=out_fi)

    
elif args["-dmr"] != "":
    out_fi =open(out_fi, "w")
    dmr = open(args["-dmr"])
    cg = [[],[],[],[],[]]
    chg = [[],[],[],[],[]]
    chh = [[],[],[],[],[]]
    whole = [[],[],[],[],[]]
    len_w, len_cg, len_chg, len_chh = 0,0,0,0
    for i in dmr.readlines()[1:]:
        i = i.strip().split("\t")
        ctx = i[-1].split(",")
        ch = int(i[1])
        st = int(i[2])
        ed = int(i[3])+1
        rng = set(range(st,ed))
        length = len(rng)
        idx = at_ch.index(ch)
        whole[idx].append(rng)
        len_w += length
        if "CG" in ctx:
            cg[idx].append(rng)
            len_cg += length
        if "CHG" in ctx:
            chg[idx].append(rng)
            len_chg += length
        if "CHH" in ctx:
            chh[idx].append(rng)
            len_chh += length
    print("Load DMR file Done")


    MakeOverlap(whole, len_w,"whole")
    MakeOverlap(cg, len_cg, "CG")
    MakeOverlap(chg, len_chg, "CHG")
    MakeOverlap(chh, len_chh, "CHH")
out_fi.close()
