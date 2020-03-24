#!/usr/bin/env python3

"""
Validate the [common format] SNP data,
output only accepted SNP data with filename plus ".2"

*From my experience, some SNP data contain multiple reference seqs (!?)
 thus this step only accepts SNP data contains single ref seq,
 & convert multi-nt symbol to separate nts.


[common format]
chr|site|ref|alt

[USAGE]
(at the folder contains multi comFmt.txt files)
> python3 a.comfmt.filter.py 


Filtering comFmt.txt SNP files and make output as comFmt.txt.2


"""

import sys, glob

fi_ls = glob.glob("*comFmt.txt")
fi_ls = sorted(fi_ls)
if len(fi_ls) == 0:
	print(__doc__)
	sys.exit(1)

for fi in fi_ls:
	name = fi
	out = name+".2"
	out = open(out,"w")
	fi = [i.strip().split() for i in open(fi)]

	accepted = 0
	multiR = 0
	multiA = 0
	corrected = 0
	notcr = 0

	acgt = ["A","T","G","C","*"]
	acgt_dic = {"Y":"C,T","R":"A,G","S":"G,C","W":"A,T","K":"T,G",
		"M":"A,C","N":"A,C,G,T","B":"C,G,T","D":"A,G,T","H":"A,C,T","V":"A,C,G"}
	for i in fi:
		ref = i[2]
		alt = i[3]
		if len(ref) == len(alt) == 1:
			if ref in acgt:
				if alt in acgt:
					accepted += 1
					print(*i,sep="\t",file=out)
				else:
					alt = acgt_dic[alt]
					if ref in alt and len(alt) == 3:
						corrected += 1
						for j in alt:
							if j != "," and j != ref:
								alt = j
								data = [i[0],i[1],ref,alt]
								print(*data,sep="\t",file=out)
					else:
						multiA += 1
			
			else:
				ref = acgt_dic[ref]
				#alt = acgt_dic[alt]
				if len(ref) > 1:
					multiR += 1

	out.close()
	print(name)
	print("#line: %s" % len(fi))
	print("accepted: %s" % accepted)
	print("multiR: %s" % multiR)
	print("multiA: %s" % multiA)
	print("corrected: %s" % corrected)
	print()
