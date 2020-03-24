#!/usr/bin/env python3

"""
This code finds SNPs in given condition,
the condition is fixed now,
while the idea will help to further coding

NOV.11th.2016

"""

import sys


#fi = "all.tag-filt.biall.vcf"
fi = sys.argv[-1]
fi = [i.strip().split() for i in open(fi)]


for i in fi:
	if i[0].startswith("##"):
		pass
	elif i[0].startswith("#"):
		header = i[9:]
		print(*header,sep=" ")
	else:
		snp = i[2]
		if snp == ".":
			snp = ":".join([i[0][-1],i[1],i[3],i[4]])
		var1 = i[9:]
		var = [i.split(":")[0] for i in var1]
		gt = [i.split(":")[1] for i in var1]
		#print(var)
		#sys.exit()
		if snp.startswith("5"):
			print(snp)
			sys.exit()
		#print(var)
		#print(gt)
		#sys.exit()
		if var[0] ==  "1/1":
			if not "1/1" in [var[1], var[2] ]:
				#if var[10] != "0/0":
				snp_anl = snp.split(":")
				if len(snp_anl[2]) == len(snp_anl[3]) == 1:
					print(snp)
"""

for i in fi:
	if i[0].startswith("#"):
		if not i[0].startswith("##"):
			ti = i[9:]
			print("BP",ti[0],ti[3],ti[5],ti[7],ti[1],ti[6],ti[4],ti[2],sep="\t")
	else:
		snp = i[2]
		info = snp.split(":")
		ch = info[0]
		bp = int(info[1])
		var = i[9:]
		var = [i.split(":")[1] for i in var]
		if ch == "1" and 2714000 <= bp <= 7628000:
			print(snp, var[0], var[3], var[5],var[7], var[1], var[6], var[4], var[2], sep="\t")
"""
