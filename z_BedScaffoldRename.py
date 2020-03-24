#!/usr/bin/env python3
"""
Scaffold name converter refering new scaffold names

since I used old type of XB genome(scaffold0001, scaffold002),
I need to change the old name to the new one (xBscaffold|bp|Rs|R3)

the reference file ("scaffold2.txt") was adopt from genome.FA file with grep ">"

15.02.12 Toshimori9 

"""



import glob

inputls = glob.glob("*bedgraph")

scfinfo = {}
for i in open("scaffold2.txt"):
	i = i.strip()
	name = i.split("|")[0][2:]
	scfinfo[name] = i


for i in inputls:
	prefix = "new/"+i
	out = open(prefix,"w")
	for j in open(i).readlines():
		j = j.strip().split("\t",1)
		scf = j[0]
		content = j[1]
		try:
			scf = scfinfo[scf]
		except:
			print(scf)
			pass
		print(scf, content, sep="\t", file=out)
	out.close()
	print("%s Done!" % i)
