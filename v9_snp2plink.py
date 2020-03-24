#!/usr/bin/env python3

"""
Generate MAP/PED file for PLINK analysis from a manipulated merged snp.txt, resultant from v7, v8_merge.comFmt_snp.py

[USAGE]
python3 snp2plink.py [option]  *snp.txt

[OUTPUT]
*snp.map
*snp.ped

[OPTION]
-pheno: a phenotype file to add the info into *snp.map file
-fill : flag this to fill the blanks with reference sequence. (default: keep blanks)


## About phenotype data, it might be better use --pheno in PLINK, with separate phenotype file (DEC10th2015)
the pheno file should be like below with tab-deliminated

=========== only header starting with "#" is acceptable (ignored in the process)
[sample_id] [pheno1] [pheno2] ....
[sample_id] [pheno1] [pheno2] ....
[sample_id] [pheno1] [pheno2] ....

Multiple *snp.ped files would be generated owing to input phenotypes number

*snp.1.ped
*snp.2.ped
...


====== Change Log ======
DEC10th 2015	Equip option to fill blanks
Nov12th 2015	Equip multi-ped file generation part
Nov02nd 2015	First build

Junesk9 

"""

import sys

##################################### Arguments

args = sys.argv[1:]
if len(args) == 0:
	print(__doc__)
	sys.exit(1)

snp_fi = args[-1]
prefix = ".".join(snp_fi.split(".")[:-1])

try:
	phe_idx = args.index("-pheno") + 1
	pheno = args[phe_idx]
	phe_mode = True
except:
	#print() #deactivated DEC10th2015
	#print("[WARNING] No phenotype data is available, the phen column (6th) in *ped file will be filled with 0")
	#print("[WARNING] You need use [--pheno] option in use of PLINK ")
	phe_mode = False
if "-fill" in args:
	fill_mode = "FILL"
else: fill_mode = "KEEP"


print("""
	
	PLINK MAP/PED generator

	INPUT     : %s
	OUT_prefix: %s .map/ped
	PHENO_DB  : %s [default: False]
	FILL_BLANK: %s [default: KEEP] (flag "-fill" to fill the blanks as ref)

""" % (snp_fi, prefix, phe_mode, fill_mode))
if phe_mode == False:
	print("[WARNING] No phenotype data is available, the phen column (6th) in *ped file will be filled with 0")
	print("[WARNING] You need to flag [--pheno] option in use of PLINK ")


snp = [i.strip().split() for i in open(snp_fi)]
samples = snp[0][5:]

################ Generate map file

map_f = prefix+".map"
map_f = open(map_f,"w")

## Since PLINK 1.9 only accept 59 chromosomes, for more than numbers, the chromosome number shoude recognized as NOT pure numeric
ch_count = []
for i in snp[1:]:
	ch_count.append(i[0])
ch_count = list(sorted(set(ch_count)))
if len(ch_count) > 59:
	exc_chr = True
	print("[WARNING] Detect total chromosome (scaffold) #: %s, excessive chromosome mode turned on" % len(ch_count))
else:
	exc_chr =  False
	print("[PROCESS] Detect total chromosome (scaffold) #: %s" % len(ch_count))


for i in snp[1:]:
	ch = i[0]
	if ch_count.index(ch) > 59: ## Masking chromosomes number more than 59, 60 convert to U60
		ch = "U"+ch
	idn = i[2]
	cM = "0"
	bp = i[1]
	print(ch,idn,cM,bp,sep="\t",file=map_f)
map_f.close()
print("[PROCESS] Generate MAP file Done")


################ Process phenotype data

phe_dic = {}
if phe_mode == False:
	for i in samples:
		phe_dic[i] = []
		phe_dic[i].append(0)
		phe_num = 1
elif phe_mode == True:
	pheno = [i.strip().split() for i in open(pheno)]
	for i in pheno:
		if not i[0].startswith("#"):
			phe_dic[i[0]] = i[1:]
			phe_num = len(i[1:])

################# Generate ped file.


for k in range(phe_num):
	file_k = str(k+1)
	if phe_mode == True:
		ped = prefix+"."+file_k+".ped"
	elif phe_mode == False:
		ped = prefix+".ped"
	ped = open(ped,"w")

	data = []
	for en, i in enumerate(samples):
		tmp = [i,i,"0","0","0"]
		if phe_mode == True:
			phe = phe_dic[i][k]
		elif phe_mode == False:
			phe = "0"
		tmp = tmp + [phe]
		data.append(tmp)
	snp = [i[3:] for i in snp[1:]]

	print("[PROCESS] Numbers of input accessions, input snp data: %s %s" % (len(samples),len(snp[0])-2))

	for i in snp:
		ref = i[0]
		i = i[2:]
		for en, j in enumerate(i):
			if j == "-":
				if fill_mode == "FILL":
					data[en].append(ref)
					data[en].append(ref)
				elif fill_mode == "KEEP":
					data[en].append("0")
					data[en].append("0") ## keep missing data as missing, expect imputation DEC9th2015
				else:
					print("[ERROR] Anadaptable option found for -fill: %s" % fill_mode)
					sys.exit(1)
			else:
				data[en].append(ref)
				data[en].append(j)

	for i in data:
		print(*i,sep="\t",file=ped)
	ped.close()
	print("[PROCESS] Generate PED file %s/%s Done" % (file_k,phe_num))
	print()
