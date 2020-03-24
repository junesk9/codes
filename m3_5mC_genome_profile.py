#!/usr/bin/env python3
"""
Survet global methylation status in view of chromosome
with a specific window size (default: 100k)
input as *input file for 5mC comparing

Three output files will emerged by the seq context (CG,CHG & CHH)

Oct. 2014 by ToshiMori
"""

import sys
import time
ctime = time.time()

########################################################
########### Args

win = 200000 # Window size for view global 5mC rate
fi = sys.stdin


##########################################################
########## Methods

def SurveyChr(data):
        ch = []
        for i in data:
                i = i.split('\t')
                ch.append(i[0])
        ch = list(sorted(set(ch)))
        return ch

def PrepareList(ch):
        CG = []
        CGout = []
        ch_max = []
        for j in ch:
                CG.append({})
                CGout.append([])
                ch_max.append(0)
        return CG, CGout, ch_max

def CalcAvg(cg, out):
        for i in ch:
                num = int(int(ch_max[ch.index(i)])/win) + 1
                for j in range(0,num):
                        out[ch.index(i)].append([0,0])
        
        for enu,k in enumerate(cg):
                for l in k.keys():
                        ind = int(l/win) # index by site info if l < 100000, ind = 0
                        out[enu][ind][0] += k[l]
                        out[enu][ind][1] += 1
        for m in out:
                for x,n in enumerate(m):
                        n[1] = str(float(n[0]/n[1]))
                        n[0] = str(x * win)
        return out


        

########################################################33


data = [i for i in fi]
ch = SurveyChr(data)
CG, CGout, ch_max = PrepareList(ch)
CHG, CHGout, ch_max = PrepareList(ch)
CHH, CHHout, ch_max = PrepareList(ch)

for k in data:
	k = k.strip().split('\t')
	loc = int(k[1])
	value = int(k[3])/(int(k[3])+int(k[4]))
	ch_max[ch.index(k[0])] = loc if loc > ch_max[ch.index(k[0])] else ch_max[ch.index(k[0])]
	if k[5] == "CG":
		CG[ch.index(k[0])][loc] = value
	if k[5] == "CHG":
		CHG[ch.index(k[0])][loc] = value
	if k[5] == "CHH":
		CHH[ch.index(k[0])][loc] = value


CGout = CalcAvg(CG, CGout)
CHGout = CalcAvg(CHG, CHGout)
CHHout = CalcAvg(CHH, CHHout)

for i in range(0,len(ch)):
                print(">Chr"+ch[i])
                for x,j in enumerate(CGout[i]):
                        try:
                                j.append(CHGout[i][x][-1])
                                j.append(CHHout[i][x][-1])
                        except: print("warning check the chr length")
                        j = str.join("\t", j)
                        print(j)
#print("Done, takes %d sec" % int(time.time() - ctime))

