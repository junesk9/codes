#!/usr/bin/env python3

"""
*******************************************************************
***** You can find more easy way by using vB_ with -ld option******
***** You dont need to run this code without specific purpose******
************Feb 22nd 2016 Junesk9**********************************
*******************************************************************

Elute SNP & R2(>=0.1) info to highlight snp-in-interest from the  manhattan plot
To gain more speed, some shell command integrated
requires 3 files in the running patha,

plink resultant ld file: generally ends with ".ld" <- indicated by command
plink resultant p-val file: generally ends with ".adjusted" <- auto-detected
snp annotation file from v4_snp_annotator, the file name ends with ".annt.csv"

[USAGE]
>cat [ld file] | grep "[snp-in-interest]" | ./plink2highlight.py [-ld|-flmm|-qqmna] -annt

[OPTION]
**One of below three is mendatory**
-ld   : PLINK p-value file (*.adjusted)
-flmm : Fast-lmm p-value file (*.txt)
-qqman: processed file for qqman() by vA_plink2qqman.py (*qqman.txt)
 

-annt : SNP annotation file by v4_plink-snp_annotator.py (*.annt.csv)
-r2   : R-square threshold, default as 0.1

[output]
[snp-in-interest].highlight.txt

====Log====
DEC17th 2015	Equip to accept various LD files (plink, fastlmm, qqman)
		Equip R-square threshold
		output include DPrime also (while it possible)
Nov4th 2015	First build

"""


import sys, glob
import math


######################### Arguments, Pre-treatment
fi = sys.stdin
""" ## DEC17th2015 Do not work, need to revise
try:
	foo = len(fi) ## Just evaluate the proper input presence
except TypeError:
	print(__doc__)
	sys.exit(1)
"""

### pre-process of STDIN
tmps = []
for i in fi:
	i = i.strip().split()
	tmp = [i[2],i[5],i[6],i[7]] ## [SNP_A,SNP_B,R2,DP]
	#print(tmp)
	#sys.exit()
	tmps.append(tmp)

inp = list(sorted(set(tmps[0]).intersection(set(tmps[1]))))[-1]
#print(list(sorted(set(tmps[0]).intersection(set(tmps[1])))))
#print(inp)
#sys.exit()
prefix = "-".join(inp.split(":"))


### Argument process
args = sys.argv[1:]
try:
	ld_idx = args.index("-ld")+1
	ld_fi = args[ld_idx]
	ld_mode = "plink"
except ValueError:
	try:
		ld_idx = args.index("-flmm")+1
		ld_fi = args[ld_idx]
		ld_mode = "flmm"
	except ValueError:
		ld_idx = args.index("-qqman")+1
		ld_fi = args[ld_idx]
		ld_mode = "qqman"
else:
	print("[ERROR] No LD file is assigned")
	print(__doc__)
	sys.exit()
	#ld_fi = glob.glob("*.adjusted")[0]

try:
	annt_idx = args.index("-annt")+1
	annt = args[annt_idx]
except:
	#annt = glob.glob("*.annt.csv")[-1]
	annt = False
try:
	r2_idx = args.index("-r2")+1
	r2_cut = float(args[r2_idx])
except:
	r2_cut = 0.1

out = prefix+".r"+str(r2_cut)+".highlight.txt"
print("""

	Highlighted SNP parser

	SNP      : %s
	LD_file  : %s (%s)
	R2_cutoff: %s (Default: 0.1)
	ANNT_file: %s
	OUTPUT   : %s
	
""" % (inp, ld_fi, ld_mode, r2_cut, annt, out))

#out = prefix+".r"+str(r2_cut)+".highlight.txt"
out = open(out,"w")
############################### Body

r2 = {inp:[inp.split(":")[1],"1","1"]} ## highlighted SNP self association data [snp,BP,R2,DP]
keys = [inp]
for i in tmps: ## [SNP_A,SNP_B,R2,DP]
	i.remove(inp)
	snp = i[0]
	bp = i[0].split(":")[1]
	r2_v = i[-2]
	dprime = i[-1]
	if float(r2_v) >= r2_cut: # r-square threshold
		tmp = [bp,r2_v,dprime]
		r2[snp] = tmp
		keys.append(snp)
print("[PROCESS] %s of SNP info will be collected" % len(keys))

print("[PROCESS] Start parse LD information from %s" % ld_fi)
for i in open(ld_fi).readlines()[1:]:
	i = i.strip().split()
	if ld_mode == "plink":
		key = i[1]
		logp = -math.log(float(i[2]),10)
	elif ld_mode == "flmm":
		key = i[0]
		try:  ## Sometimes Fast-lmm outputs p-value as 0 or "NaN"
			logp = -math.log(float(i[4]),10)
		except ValueError:
			if i[4] == "0":
				logp = 20
			elif i[4] != 0:
				logp = "NA"
	elif ld_mode == "qqman":
		key = i[0]
		try:  ## Sometimes Fast-lmm outputs p-value as 0 or "NaN"
			logp = -math.log(float(i[3]),10)
		except ValueError:
			if i[3] == "0":
				logp = 20
			elif i[3] != 0:
				logp = "NA"		
	if key in keys: 
		r2[key].append(logp)
## Validation stage
non_p = 0
for key in r2.keys():
	if len(r2[key]) != 4: #[BP,R2,DP,logP]
		if ld_mode in ["qqman","flmm"] and len(r2[key]) == 3:
			r2[key].append("NA") # Assign the missing pval as NA
			#print("[WARNING] %s do not have any P-val, assigned as 0" % key) 
			non_p += 1
		else:
			print("[ERROR] unknone error in parse logP value")
			sys.exit(2)
print("[WARNING] %s SNPs do not have any P-value, assigned as 0" % non_p)

if annt != False:
	print("[PROCESS] Start parse SNP annotations from %s" % annt)
	annt = [i.strip().split(",") for i in open(annt)]
	for i in annt:
		if i[1] in keys:
			desc = i[-1]
			key = i[1]
			r2[key].append(desc)
	print("SNP","BP","R2","DP","logP","ANNT",sep="\t",file=out)
else: # annt == False
	print("[PROCESS] Skip the SNP annotation step, no ANNT file assigned")
	print("SNP","BP","R2","DP","logP",sep="\t",file=out)


### Make output file
print("[PROCESS] Make the output file")
for key in sorted(keys):
	print(key,end="\t",file=out)
	print(*r2[key],sep="\t", file=out)
out.close()


