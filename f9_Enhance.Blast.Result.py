#!/usr/bin/env python3

"""
Add hit_len/query_len to the local blast output file (format6)

mandatory options are below
-b: the blast output file to add info
-q: the query multi-fasta file

the output file will be generated as ends with "enhanced.txt"


==log==
Jan 6th 2016	First build


Junesk9  
"""

import sys

args = {"blast":"","query":""}

arg = sys.argv[1:]
try:
	blast_idx = arg.index("-b")+1
	query_idx = arg.index("-q")+1
except:
	print(__doc__)
	sys.exit(2)
args["blast"] = arg[blast_idx]
args["query"] = arg[query_idx]

blast = args["blast"]
query = args["query"]
prefix = ".".join(blast.split(".")[:-1])
out = prefix+".enhanced.txt"

###################################################

q = {}
for i in open(query):
	i = i.strip()
	if i.startswith(">"):
		name = i.split()[0][1:]
		q[name] = []
	else:
		q[name].append(i)

for key in q.keys():
	q[key] = "".join(q[key])

#################################################

out = open(out,"w")
header = ["query","hit","identity","hit_len","mismatch","gap","query_st","query_ed","hit_st","hit_ed","E-val","score","hit/query"]
print(*header,sep="\t",file=out)
for i in open(blast):
	i = i.strip().split()
	query = i[0]
	query_len = len(q[query])
	hit_len = int(i[3])
	query_oc = round(hit_len/query_len*100,2)
	i.append(query_oc)
	print(*i,sep="\t",file=out)
out.close()
