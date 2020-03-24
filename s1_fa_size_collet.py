#!/usr/bin/env python3

"""
Collect the reads with correct size from a multi-fasta file
the multiple read sizes can be assigned, and out as separate files

USAGE:
python3 s1_fa_size_collet.py [.fa]


Jan. 22nd. 2015 ToshiMori 

"""

import sys, time

### Argumnts
inp = sys.argv[-1]
inp = [i.strip() for i in open(inp)]
prefix = inp.split(".")[0]
wantedSize = [21,22,24]

### Body
ctime = int(time.time())
dic = {}
for en,i in enumerate(inp):
	if en % 2 == 0:
		i = i.split(" ")[0]
		key = i
		dic[key] = ""
	elif en % 2 == 1:
		value = i
		dic[key] = value
rtime = int(time.time()) - ctime
print("load %s done! takes %s sec" % (inp,rtime))

### Function
def SizeCollect(dic, size):
	out_name = prefix+".size_"+str(size)+".fa"
	out = open(out_name, "w")

	k = 0
	for i,j in dic.items():
		if len(j) == int(size):
			print(i,file=out)
			print(j,file=out)
			k += 1

	rtime = int(time.time()) - ctime
	print("make %s file Done! takes %s secs" % (out_name, rtime))
	out.close()
#end of def

for x in wantedSize:
	SizeCollect(dic,x)
