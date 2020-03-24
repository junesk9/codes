#!/usr/bin/env python3

"""
Restriction enzyme cut fragment extractor 
for RAD-seq analysis

[USAGE]
python3 vB...py -dedup  [fa]

Junesk9

==Change log===
Jun 02nd 2016	Equip de-duplication step, the same process to fC_dedup.fa.
Jun 01st 2016	Bug fixed, applicapable to differ-length RE set
May 31th 2016	First Build

"""

import sys 
import time
ctime = time.time()

args = sys.argv[1:]
if len(args) == 0:
	print(__doc__)
	sys.exit()
dedup = True if "-dedup" in args else False

RE = {"AGATCT": "BglII",
      "GAATTC": "EcoRI",
      "CTGCAG": "PstI",
      "CCGG": "HpaII"}


fi = sys.argv[-1]
re2 = "CCGG" #PstI #"AGATCT" ## BglII
re1 = "CTGCAG" #HpaII  #"GAATTC" ## EcoRI
min_s = 100
max_s = 500
prefix = "wheat.genome.RE-cut."+str(int(ctime))+".fa"
dedup = True if "-dedup" in args else False

print("""

	====== RE-cut genome preparation ======

	input      : %s
	RE1        : %s [%s]
	RE2        : %s [%s]
	retain_rng : %s - %s nt
	output     : %s
	
	dedup_mode : %s 

""" % (args[-1], re1, RE[re1], re2, RE[re2], min_s, max_s, prefix, dedup))


genome = {}
for i in open(fi):
	i = i.split()[0]
	if i.startswith(">"):
		name = i[1:]
		genome[name] = []
	else:
		genome[name].append(i)

for ch in genome.keys():
	genome[ch] = "".join(genome[ch])

print("[PROGRESS] Genome loading done, takes %s min" % int((time.time()-ctime)/60))
print("[PROGRESS] %s chromosomes/scaffolds are detected" % len(genome.keys()))

out = open(prefix,"w")
x, y, z = 0, 0, 0
for ch in genome.keys():
	x += 1
	re1_site, re2_site = [],[]
	for en,nt in enumerate(genome[ch]):
		if en+6 < len(genome[ch]):
			re1_len = len(re1)
			re2_len = len(re2)
			ctx1 = genome[ch][en:en+re1_len]
			ctx2 = genome[ch][en:en+re2_len]			
			ctx1, ctx2 = ctx1.upper(), ctx2.upper()
			if ctx1 == re1:
				re1_site.append(en)
			elif ctx2 == re2:
				re2_site.append(en)
	for re1s in re1_site:
		k = max_s + 1
		border = []
		for re2s in re2_site:
			dist = abs(re1s - re2s)
			if min_s <= dist < k :
				k = dist
				border = [re1s, re2s+re2_len] if re1s < re2s else [re2s, re1s+re1_len]
		#end of for
		if not border == []:
			head = ch.split()[0].split("_")[2:][::-1] ## extract from TGACv1_scaffold_018948_1AL to [1AL, "018948"]
			#print(head, border)
			head = ">"+"_".join(head)+"_"+str(border[0])+":"+str(border[1]) ## >1AL_018948_st:ed
			seq = genome[ch][border[0]:border[1]]
			print(head,seq,sep="\n", file=out)
			z += 1 
		#end of if			
	if x > 9999:
		x = 0
		y += 10000
		print("[PROGRESS] %s chrs/scaffolds processed, %s seqs retained, takes %s min" % (y, z,int((time.time()-ctime)/60)))
	else: pass
	#end of if
#end of for
print("[PROGRESS] RE-cut fragment extraction done, %s sequences have been retrived totally taked %s min" % (z, int((time.time()-ctime)/60)))
out.close()

def revcomp(fa):
        tmp = ""
        for f in fa[::-1]:
                if f.upper() == "A":
                        tmp += "T"
                elif f.upper() == "T":
                        tmp += "A"
                elif f.upper() == "G":
                        tmp += "C"
                elif f.upper() == "C":
                        tmp += "G"
                else:
                        tmp += f
        return tmp


if dedup == True:
	from collections import Counter
	fi = prefix
	print("[PROGRESS] Start de-duplication")
	ctime = time.time()

	dic, dic_r, dic2 = {}, {}, {}
	x, n = 0, 0
	for i in open(fi):
		i = i.strip()
		if i.startswith(">"):
			name = i
			dic[name] = ""
			dic_r[name] = ""
			x += 1
		else:
			dic[name] = i.upper()
			dic_r[name] = revcomp(i)
			dic2[dic[name]] = name
		if x > 99999:
			x = 0
			n += 100000
			print("[PROGRESS] %s seqs are loaded, takes %s sec" % (n, int(time.time()-ctime)))
	print("[PROGRESS] loading fasta done, %s seqs are loaded, takes %s min" % (n+x, int((time.time()-ctime)/60)))

	if len(dic.keys())*2 == len(set(list(dic.values())+list(dic_r.values()))):
		print("[DONE] No duplicate is found; All progress done")
	else:
		print("[PROGRESS] Duplication was detected, de-duplication step activated")
		seqs = list(dic.values()) + list(dic_r.values())
		count = Counter(seqs)
		y, z, m = 0,0,0
		prefix2 = ".".join(fi.split(".")[:-1])+".dedup.fa"
		out2 = open(prefix2, "w")
		for i,j in count.items():
			y += 1
			if j == 1:
				try:
					name = dic2[i]
					print(name, i,sep="\n", file=out2)
					m += 1
				except KeyError: pass
			else: pass
			if y > 99999:
				y = 0
				z += 100000
				print("[PROGRESS] %s seqs are processed, takes %s min" % (z, int((time.time()-ctime)/60)))
			else: pass
		out2.close()
		seqs, dic, dic_r, dic2 = [],[],[],[]
		print("[DONE] Unique %s seqs retained from %s seqs" % (m, n+x))


