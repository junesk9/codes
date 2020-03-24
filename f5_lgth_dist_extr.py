#!/usr/bin/env python3

"""
Multi-fasta length distribution extractor

USAGE:
python3 histo.py [fa1] (optional [fa2] [fa3] ...)
(if multiple files loaded, the merged data will be out.)

OUTPUT:
STDOUT
[len][fraction][count]; Tab-delimited

Aug. 23rd 2015 Junesk9

"""

import sys

#################### Arguments

args = sys.argv[1:]
if len(args) == 0:
	print("""
################# Multi-fa length distribution Extractor ######

Usage: 
python3 histogram.py [fa1] (optional [fa2] [fa3] ...)

""")
	sys.exit()

bin_size = 100
if "-bin" in args:
	idx = args.index("-bin")
	bin_size = args[idx+1]
	if str(bin_size).isnumeric() == False:
		print("[ERROR] The bin size should be integer!")
		sys.exit()

################### Body

dic = {}
for fi in args:
	fi = open(fi)
	for i in fi:
		i = i.strip()
		if i.startswith(">"):
			name = i
			dic[name] = []
		else:
			dic[name].append(i)

lgth = []
for i in dic.keys():
	dic[i] = "".join(dic[i])
	lgth.append(len(dic[i]))


rng = [i for i in range(0,max(lgth)+bin_size,bin_size)]
hist = []
for j in rng:
	hist.append(0)

for i in lgth:
	idx = int(round(i,-2)/bin_size) ## distill the index at "hist" list eg) 11345 -> 11 (the 12nd location of "hist")
	hist[idx] += 1 # count-up at the "hist" list

tot = sum(hist)
tot_pr = 0
print("lgth","portion","count",sep="\t")
for en,i in enumerate(rng):
	if i != 0:
		tot_pr += float(hist[en]/tot)
		print(i, round(float(hist[en]/tot),3) ,hist[en], sep="\t")
print("total",str(round(tot_pr,2)),str(tot),sep="\t")
	
		
