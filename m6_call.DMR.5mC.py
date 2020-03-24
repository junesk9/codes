#!/usr/bin/env python3
"""
Call the average 5mC level by seq. context of each allocated region (generally DMRs) 
from the genmoe 5mC files (*input format).

[Usage]
python3 m6_call.DMR.5mC.py [region.file]

*the multiple genome 5mC file inputs are assigned by folder (glob)
*[region.file] must contain 5-columns with allocated info (common dmr output format). 
[region.name][chr][st][ed] ........   [context]

Nov.2014 Toshimori

[CHANGE LOG]
MAY 3rd 2015	Run time expression changes: secs -> ca. hrs

MAY 2nd 2015	Also acceptable the intact *gff file as [region.file]
                equip the conversion part
"""
import time, os, sys, pickle, glob, shutil
from multiprocessing import Process


ctime = time.time()
fi_ls = glob.glob('../*input')
dmr = sys.argv[1]
threads = 3 # Set the maximum No. of threads to be used

##########################################################
######### Methods

def find_5mC(dmr, fi, count, ctx):
    ls = [[],[],[]]
    for j in dmr:
        for k in fi:
            if int(j[1]) <= int(k[1]) <= int(j[2]):
                if k[3] == 'CG':
                    ls[0].append(k[2])
                elif k[3] =='CHG':
                    ls[1].append(k[2])
                elif k[3] =='CHH':
                    ls[2].append(k[2])
    pick = folder+"/"+count+"."+ctx+".out.pickle"
    with open(pick, 'wb') as f:
        pickle.dump(ls, f)

def multiproc(ch, n, threads, fi_name):
    if n == 0:
        na = "CG"
    elif n == 1:
        na = "CHG"
    elif n == 2:
        na = "CHH"

    ## Optional parts: If too many (No.>threads) chromosomes (scaffolds) should be handled
    ch2 = []
    x, y = 0, threads
    if len(ch) > threads:
        print("Warning!, The input %s harbors %s of scaffolds(chromosomes), while %s of threads is set to use" % (fi_name, len(ch), threads))
        while len(ch) > x:
            ch2.append(ch[x:y])
            x += threads
            y += threads
    else: ch2 = [ch]
    #### mar.16.15
    
    for en1, chrs in enumerate(ch2):
        p_ls = []
        for en2, i in enumerate(chrs):
            ch_idx = str(en1)+str(en2)
            pick = folder+"/"+i+str(n)+".input.pickle"
            pick2 = folder+"/"+i+str(n)+".dmr.pickle"
            with open(pick, 'rb') as f:
                fi = pickle.load(f)
            with open(pick2, 'rb') as f:
                dmr = pickle.load(f)
            p = Process(target=find_5mC, args=(dmr, fi, ch_idx, na))
            p.start()
            p_ls.append(p)
                  
        for p in p_ls:
            p.join()
            rtime = int(time.time() - ctime)
            if rtime > 3600: ## Convert time expression secs -> hrs
                rtime = str(round(rtime/3600, 2)) + " hrs"
            else:
                rtime = str(rtime) + " secs"
            print("process done takes %s: %s %s" % (rtime, fi_name, na))



################################################################################
############### Body

ch = []
dmr_ls =[[],[],[]]

'''
##Convert .gff to .dmr
if len(dmr[0]) > 5:
    print("Warning!! You inputed a GFF3 file than the 5mC.input file!")
    tmp = ["void-not-meaning"]
    for en, i in enumerate(dmr):
        name = str(en)
        chrs = i[0].split("|")[0][2:]
        st = i[3]
        ed = i[4]
        ctx = "CG,CHG,CHH"
        tmp.append([name,chrs,st,ed,ctx])
    dmr = tmp
    tmp = []
else: pass
## may.02.15
'''

for i in open(dmr).readlines()[1:]:
    i = i.rstrip().split('\t') 
    ch.append(i[1])
    if "CG" in i[-1]:
        dmr_ls[0].append(i[1:4])
    if "CHG" in i[-1]:
        dmr_ls[1].append(i[1:4])
    if "CHH" in i[-1]:
        dmr_ls[2].append(i[1:4])
ch = sorted(set(ch))

for fi in fi_ls:
    folder = 'tmp_'+str(int(time.time()))
    os.makedirs(folder)
    #print(fi, dmr)
    #sys.exit()
    name = fi.split("/")[-1].split(".")[0]+"_"+dmr.split("/")[-1].split(".")[1]+".5mC.txt"
    #name = fi.split("/")[-1].split(".")[0]+"_"+str.join('.',dmr.split(".")[:-2])+".5mC.txt"
    out = open(name, 'w')

    ## Parse the .input data by seq context
    fi_ctx = [[],[],[]]
    for i in open(fi).readlines():
        i = i.split('\t')
        ctx = i[5]
        
        try: ## some 5mC input linse contain casual context info, especially from bwa-meth
            rate = float(int(i[3])/(int(i[3])+int(i[4])))
            if ctx == "CG":
                fi_ctx[0].append([i[0],i[1],rate,ctx])
            elif ctx == "CHG":
                fi_ctx[1].append([i[0],i[1],rate,ctx])
            elif ctx == "CHH":
                fi_ctx[2].append([i[0],i[1],rate,ctx])
        except ValueError: pass
    
    
    for n, d in enumerate(dmr_ls):
        # Parse DMR input data to list [tmp], by chromosome
        tmp = []
        for i in range(len(ch)):
            tmp.append([])

        for i in d:
            idx = ch.index(i[0])
            tmp[idx].append(i)

        for en, i in enumerate(ch):
            pick = folder+"/"+i+str(n)+".dmr.pickle"
            with open(pick, 'wb') as f:
                pickle.dump(tmp[en], f)

        # Parse 5mC input data to list [tmp], by chromosome
        tmp = []
        for i in range(len(ch)):
            tmp.append([])

        for i in fi_ctx[n]: # Some input files contains other chr No. (6,7..)
            try:
                idx = ch.index(i[0])
                tmp[idx].append(i)
            except ValueError: pass

        for en, i in enumerate(ch):
            pick = folder+"/"+i+str(n)+".input.pickle"
            with open(pick, 'wb') as f:
                pickle.dump(tmp[en], f)
        tmp = []  # blank the tmp
        fi_name= fi # For enhanced verbose Apr.6th,2015
        multiproc(ch, n, threads, fi_name)
    #end of for

    out_ls =[[],[],[]]
    pick_ls = glob.glob(folder+"/*out.pickle")
    for i in pick_ls:
        with open(i, 'rb') as f:
            data = pickle.load(f)
            out_ls[0] += data[0]
            out_ls[1] += data[1]
            out_ls[2] += data[2]

    print("CG","CHG","CHH",sep="\t",file=out)
    max_num = max(len(out_ls[0]), len(out_ls[1]), len(out_ls[2]))
    for l in range(0,max_num):
        try: print(out_ls[0][l], end="\t", file=out)
        except: print('', end="\t", file=out)
        try: print(out_ls[1][l], end="\t", file=out)
        except: print('', end="\t", file=out)
        try: print(out_ls[2][l], file=out)
        except: print('', file=out)
    #end of for
    out.close()
    shutil.rmtree(folder)

