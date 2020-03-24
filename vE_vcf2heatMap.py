#!/usr/bin/env python3

import sys


vcf = sys.argv[-1]
prefix = ".".join(vcf.split(".")[:-1])
out = prefix +".heatmap.txt"

out = open(out,"w")
for i in open(vcf):
	if i.startswith("##"):
		pass
	elif i.startswith("#C"):
		i = i.strip().split()
		idx = i[:3]+i[9:]
		print(*idx,sep="\t",file=out)
	else:
		i = i.strip().split()
		idx = i[:3]
		snp = []
		for j in i[9:]:
			j = j.split(":")[0].split("/")
			try:
				j = int(j[0])+int(j[1])
			except ValueError:
				j = "NA"
			snp.append(j)
		idx = idx + snp
		print(*idx,sep="\t",file=out)
out.close()

		
