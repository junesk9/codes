#!/usr/bin/env python3

"""

slice fa with chr:st-ed set range

tmp.fa-slice.py tair10.fa 1:10000-20000 [2:100000-300000]]

DEC.08 2016 Junesk9
"""


import sys


args = sys.argv[1:]
fa = args[0]
rngs = args[1:]
for r in rngs:
	if not ":" in r:
		print("[ERROR] the range would be given as chr:st-ed forms")
		sys.exit()

dic = {}
for f in open(fa):
	f = f.strip()
	if f.startswith(">"):
		name = f.split()[0][1:]
		dic[name] = []
	else:
		dic[name].append(f)
chrs = sorted(list(dic.keys()))
for k in chrs:
	dic[k] = "".join(dic[k])
fa = dic


for r in rngs:
	ch = r.split(":")[0]
	st = int(r.split(":")[1].split("-")[0])
	ed = int(r.split("-")[-1])
	if not ch in chrs:
		print("[ERROR] Cannot find chromosome names: %s" %  ch)
		print("[ERROR] Only possible chrmosomes are: %s" % ",".join(chrs))
	else:
		name = ">"+r
		seq = fa[ch][st-1:ed]
		print(name)
		for i in range(0,len(seq),60):
			print(seq[i:i+60])


