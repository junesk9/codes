#!/usr/bin/env python3

import sys

######################argument
args = sys.argv
win = 1000000

##################parsing genome data
arg = {"g":"", "gff":"", "i":""}
if "-g" in args:
    gidx = args.index["-g"]+1
    arg["g"] = args[gidx]

    gdic = {}
    gn = arg["g"]
    for i in open(gn):
        i = i.strip()
        if i.startswith(">"):
            ch = i.split()[0][1:]
            gdic[ch] = []
        else:
            gdic[ch] + [i]
    for i in gdic.keys():
        gdic[i] = "".join(gdic[i])
        gdic[i] = len(gdic[i])
else:
    gdic = {"2":19698289,
            "1": 30427671,
            "4": 18585056,
            "5": 26975502,
            "3": 23459830} #TAIR10 genome 

hist = {}
for k in gdic.keys():
    hist[k] = []
    win_size = 1+int(int(gdic[k])/win)
    for x in range(win_size+1):
        hist[k].append(0)

##################### parsing gff
 
if "-gff" in args:
    gidx = args.index["-gff"]+1
    arg["gff"] = args[gidx]
else:
    arg["gff"] = "tair10_m.gff"

gff = arg["gff"]
tmp = {}
whole =False
for g in open(gff):
    if g.startswith("#"):
        pass
    else:
        g = g.strip().split()
        cls = g[2]
        if cls == "mRNA" :
            agi = g[-1].split(".")[0].split("=")[-1]
            if agi in tmp.keys():
                pass
            else:
                ch = g[0]
                st = g[3]
                tmp[agi] = [ch, st]
            #end of if
            if whole == True:
                st = int(int(st)/win)
                try:
                    hist[ch][st] += 1
                except KeyError:
                    pass
        else: pass
        #end of if
    #end of if
#end of for
gff = tmp

if whole == True:
    keys = sorted(list(hist.keys()))
    for k in keys: 
        print(k,end="\t")
        print(*hist[k],sep="\t")
    sys.exit()
    
##################### distribution of given gene list
if "-inp" in args:
    gidx = args.index["-gff"]+1
    arg["i"] = args[gidx]
else:
    arg["i"] = "input.txt"
 
inp = arg["i"]
inp = [i.strip().split(".") [0] for i in open(inp)]
inp = list(set(inp))

x, y = 0, 0
for i in inp:
    y += 1
    try:
        ich = gff[i][0]
        ist = gff[i][1]
        ist = int(int(ist)/win)
        hist[ich][ist] += 1
    except KeyError:
        x += 1
        

keys = sorted(list(hist.keys()))
print("[RESULT] Totally %s attemped, %s  rejected (e.g. transposon)" % (y,x))
for k in keys:
    print(k,end="\t")
    print(*hist[k],sep="\t") 

    

