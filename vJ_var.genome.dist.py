#!/usr/bin/env python3
"""
Genereate Chromosomal distribution of VCF file

python [-g: genome.fasta] VCF


"""

import sys

######################argument
args = sys.argv
win = 500000
th = 0.9 ## Threshold for counting individual data of merged VCF file
ems = True

##################parsing genome data
def makeHist():
	arg = {"g":""}
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

	return hist
###############parsing vcf

vcf = args[-1]
ch, info = [], []
vcf = [i.strip() for i in open(vcf)]
for i in vcf:
	if not i.startswith("##"):
		if i.startswith("#CHROM"):
			ind = i.split()[9:]
		else:
			i = i.split()
			ch.append(i[0])
			info.append(i[:2])
#print(ind)
#sys.exit()

#################Warning parts via compare chromosome indicators
hist = makeHist()
ch = set(ch)
keys = set(hist.keys())
cmn = ch.intersection(keys)

if len(cmn) == 0:
	print("[ERROR] No chromosome name matched between VCF and GENOME info")
	sys.exit()
#print(len(cmn),len(ch),len(keys))
if len(cmn) < len(ch):
	ch = ",".join(sorted(list(ch)))
	print("[WARNING] Only %s CHR data can be parsed, yout input CHR is %s" % (len(cmn), ch))


##############body


for i in info:
	ch = i[0]
	nt = int(i[1])
	if ch in hist.keys():
		idx = int(nt/win)
		hist[ch][idx] += 1

keys = sorted(list(hist.keys()))
print("Total Variants")
for k in keys:
	print(k,end="\t")
	print(*hist[k],sep="\t")


for en,i in enumerate(ind):
	hist = makeHist()
	en = en+9
	print()
	print(i)
	for v in vcf:
		if not v.startswith("#"):
			v = v.split()
			ct = v[3]+":"+v[4]
			ct = True if ct in ["C:T","G:A"] else False 
			info = v[en].split(":")[1].split(",")
			info = float(int(info[1])/(int(info[0])+int(info[1]))) if info != ['0','0'] else 0
			ler = v[-1].split(":")[1].split(",")
			ler = float(int(ler[1])/(int(ler[0])+int(ler[1]))) if ler != ['0','0'] else 0
			if info >= th:
				ch = v[0]
				nt = int(v[1])
				if ct == True or ems == False:
					if ler > 0.9:
						pass
					else:
						idx = int(nt/win)
						hist[ch][idx] += 1
				elif ems == True and ct == False:
					pass
	for k in keys:
		print(k,end="\t")
		print(*hist[k],sep="\t")


