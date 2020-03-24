#!/usr/bin/env python3

"""
Convert & concactermerize three .wig files (seperated by seq context) from BSSEQ analysis
to single bed file for further analyses

== Change log
14.11.19 - Add [ctx] verifying body, some files use mCG, mCHG, mCHH than CG, CHG, CHH
14.11.01 - Automatically detect the float length, and apply it to proper Flaction (see var: rt_ln)
14.10.30 - Change the [ClacFlac] from 'float()' to 'round()'
"""


import sys
import glob
import time

ctime = time.time()
ls = glob.glob("GSM1419610*")

ctx = [i.split('_')[-1].split('.')[0] for i in ls]
name = [i.split('_')[1] for i in ls]
fi_name = name[0]+".5mC.input"
out = open(fi_name, 'w')

# Some files annotate its sequence contetx as mCG, mCHG, mCHH than CG, CHG, CHH
for en,i in enumerate(ctx):
	if i.startswith("m"):
		ctx[en] = i[1:]

def CalcFlac(rate, rt_ln):
	mC, aC = 0, 4
	while rate != round(mC/aC, rt_ln):
		mC += 1
		if mC > aC:
			mC = 1
			aC += 1
	else: pass
	NmC = aC - mC
	return mC, NmC

k = 0
for i in ls:
	for j in open(i).readlines():
		j = j.strip()
		if j.startswith('vari'):
			chro = int(j[-1])
		elif not j.startswith('track') and j != '':
			j =  j.split("\t")
			if float(j[1]) < 0:
				st = '-'
			elif float(j[1]) >0:
				st = '+'
			elif float(j[1]) == 0:
				st = 'n'
			rate = abs(float(j[1]))
			if rate != 0:
				rt_ln = len(str(rate))-2
			else: rt_ln = 0
			site = j[0]
			mC, NmC = CalcFlac(rate, rt_ln) 
			print(chro, site, st, int(mC), int(NmC), ctx[k], "NNN", sep="\t", file=out)
	k += 1
out.close()
print("Done, takes %d secs" % round(time.time()-ctime))

