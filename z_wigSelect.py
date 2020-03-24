#!usr/bin/env python3
"""
Collect desired scaffold wig info from a whole big file
to facilitate local genome browser work

[USAGE]
python3 wigSelect.py [wig-prefix]

Jan.3rd.2015 Toshimori
"""

import glob, sys


fi = glob.glob(sys.argv[-1])
scf = ["scaffold1", "scaffold15", "scaffold45", "scaffold72", "scaffold9", "scaffold96"]

k = 1
for i in fi:
	wi_dic = {}
	prefix = ".".join(i.split(".")[:-2])+".select.wig"
	out = (prefix, "w")
	wig = open(i)
	for n,j in enumerate(wig):
		j = j.strip()
		if n == 0:
			print(j, file=out)
		else:
			if j.startswith("variable"):
				ch = j.split("=")[-1]
				if ch in scf:
					name = j
					wi_dic[name] = []
					k = 2
				else: k = 1
			elif k == 2:
				wi_dic[name].append(j)

	for i,j in wi_dic.items():
		print(i, file=out)
		for k in j:
			print(k, file=out)
	out.close()

