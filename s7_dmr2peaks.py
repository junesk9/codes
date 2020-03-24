#!usr/bin/env python3
"""
Find DMRs harbouring sRNA peaks

USAGE:
python3 dmr2peaks.py [dmr] [peaks] > STDOUT

Feb.2015 ToshiMori
"""

import sys

args = sys.argv[1:]
dmrFi = args[0]
peakFi = args[1]

dmr =[]
for i in open(dmrFi).readlines()[1:]:
	i = i.strip().split("\t")
	name = i[0]
	ch = i[1]
	st = i[2]
	ed = i[3]
	ln = i[4]
	ctx = i[-1]
	dmr.append([name,ch,st,ed,ln,ctx])

peak = []
for j in open(peakFi):
	j = j.strip().split("\t")
	name = j[0]
	ch = j[1]
	st = j[2]
	ed = j[3]
	ln = j[4]
	peak.append([name,ch,st,ed,ln])

print("#dmr","dmr_chr","dmr_st","dmr_ed","peak","peak_st","peak_ed","comm_len","dmr_ctx",sep="\t")
for x in dmr:
	for y in peak:
		if x[1] == y[1] or x[1]+"|" in y[1] or y[1]+"|" in x[1]:
			rng1 = {i for i in range(int(x[2]),int(x[3])+1)}
			rng2 = {i for i in range(int(y[2]),int(y[3])+1)}
			comm = sorted(list(rng1.intersection(rng2)))
			if len(comm) != 0:
				print(x[0],x[1],x[2],x[3],y[0],y[2],y[3],len(comm),x[-1],sep="\t")
			else: pass
		else: pass
