#!/usr/bin/env python3

"""
Conditional SNP selection for EMS F2 bulked sequence

OkaL is from Ler seq (PMID: 27354520)
Other Oka[1..15] if EMS, Ler-crossed, F2 bulk (25 ind?)

Condition is simple:
SNP ratio of EMS should be bigger than 0.8 (80%) and OkaL
SNP ratio has significant difference to OkaL (Fisher 2-exact test)
or SNP ratio 1.0 (100%) retained.


Feb 15th 2017 Junesk9


############ Log
18 03 12 Add options; -annt: if vcf.annt file given, the output contains SNP annotation column
                      -rt : the threshold of SNP ratio considered to be anaylzed (default: 0.9)
                      -neg : Indicate the negative snp data column, if not given, only 100% SNPs outputs
"""


import sys
import pandas as pd
import time
import scipy.stats as stats


################## Options

snp = sys.argv[-1]

neg = "1_NS"
th_rt = 95 ## Percentage
annt = ""
annt_mode = False
out = ".".join(snp.split(".")[:-1])+".EMS-chosen.txt"


args = sys.argv[1:]
if "-annt" in args:
	an_idx = args.index("-annt") + 1
	annt = args[an_idx]
	annt_mode = True
if "-rt" in args:
	rt_idx = args.index("-rt") + 1
	th_rt = int(args[rt_idx])
if "-neg" in args:
	neg_idx = args.index("-neg") + 1
	neg = args[neg_idx]
if "-vcf" in args:
	vcf_idx = args.index("-vcf") + 1
	snp = args[vcf_idx]
 



################## Body

#### if annt_mode == True

if annt_mode == True:
	print("[PROGRESS] SNP ANNOTATION IS GIVEN: %s" % annt)
	tmp = {}
	for i in open(annt).readlines()[1:]:
		i = i.strip().split("\t")
		tmp[i[1]] = i[2]
	annt = tmp
else:
	print("[PROGRESS] No SNP ANNOTATION IS AVAILABLE")

####
snp = [i.strip() for i in open(snp)]
for s in snp:
	if s.startswith("##"):
		pass
	elif s.startswith("#CHROM"):
		header = s.split()[9:]
		vcf = {h:[] for h in header}
		vcf["name"] = []
	else:
		s = s.split()
		name = s[2]
		vcf["name"] += [name]
		s = s[9:]
		for en,i in enumerate(s):
			vcf[header[en]] += [i]
print("[PROGRESS] Loading VCF file done")

def emsperc(string):
	## string would be int,int form
	s = string.split(",")
	s[0], s[1] = int(s[0]),int(s[1])
	try:
		perc = round((s[1]*100)/(s[0]+s[1]),2)
		return perc
	except ZeroDivisionError:
		if string == "0,0":
			return "NA"
		else:
			print("[ERROR] Zero Division Error %s" % string)
			sys.exit()
def fisher(str1, str2):
	## string would be int,int form
	str1 = str1.split(",")
	str2 = str2.split(",")
	str1 = [int(i) for i in str1]
	str2 = [int(i) for i in str2]
	
	p = stats.fisher_exact([str1,str2], alternative='less')[-1]
	return p


keys = sorted(list(vcf.keys()))
keys.remove(neg)
keys.remove("name")
print(keys)
#sys.exit()

#out = "Oka-EMS.2nd.txt"
out = open(out,"w")

header = ["GERM","SNP_ID","EMS%","NEG%","PVAL"]
if annt_mode == True:
	header.append("ANNT")
else: pass
print(*header, sep="\t",file=out)
#columns = ["GERM","SNP_ID","EMS%","LER%","PVAL"]
#df = pd.DataFrame(columns=columns)
for k in keys:
	snps = vcf[k]
	z, y = 0, 0 ## SNP counter
	#print(k, snps[:5])
	for en,snp in enumerate(snps):
		name = vcf["name"][en]
		#print(k, name, snp)
		#sys.exit()
		snp = snp.split(":")[:2]
		snpL = vcf[neg][en].split(":")[:2]
		emsp = emsperc(snp[1])
		lerp = emsperc(snpL[1])
		p = "NA" ## Temperary value
		output = [k, name, emsp, lerp, p]
		y += 1
		if annt_mode == True:
			try:
				desc = annt[name]
			except KeyError:
				print("[WARN] No annotation for SNP %s was found" % name)
				desc = "-"
			output.append(desc)

		if not "NA" in [emsp, lerp]:
			p = fisher(snp[1],snpL[1])
			output[-2] = p 
			if p < 0.05 and emsp >= th_rt > lerp:
				z += 1
				print(*output, sep="\t", file=out)
		elif emsp == 100:
			z += 1
			print(*output, sep="\t", file=out)
		"""
		if snp[0] == "1/1" and snpL[0] != "1/1":
			emsp = emsperc(snp[1])
			lerp = emsperc(snpL[1])
			print(k,name,emsp,lerp,sep="\t",file=out)

		"""
	print("[PROGRESS] Processing %s SNPs done, %s/%s SNPs eluted" % (k, z, y))
out.close()


"""
stime = time.time()
x, y = 0,0
for s in snp:
	if s.startswith("##"):
		pass
	elif s.startswith("#CHROM"):
		header = s.split()[9:]
		#print(*header, sep="\t")
		df = pd.DataFrame(columns=header)
	else:
		line = s.split()[9:]
		idx = [s.split()[2]]
		#print(idx)
		line = pd.DataFrame([line], index=idx, columns=header)
		df = df.append(line)
		x += 1
		if x >= 10000:
			y += 1
			x = 0
			ctime = int(time.time()-stime)
			print("[PROGRESS] PANDAS read %s lines, takes %s sec" % (y*10000,ctime))
print(df)
"""
