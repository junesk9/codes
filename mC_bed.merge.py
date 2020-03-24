#!/usr/bin/env python3
"""
merging two bed files from bwa-meth pathway
(Because bwa-meth occurs many errors in process with large seq files;
 Sometimes need to separate seq data)

Use Dictionary func. of python, will takes several hours

Many [counting body] equipped, since (I suspect) >20 min without output,
python3 go being mocked.

USAGE: python3 bed.merge.py [bed1] [bed2]

25th Nov. 2014 ToshiMori
"""


import sys, time

ctime = time.time()
args = sys.argv[1:]
fi1 = open(args[0])
fi2 = open(args[1])
name = args[0].split(".")[0]+".merged.bed"
out = open(name,'w')
header = "#chrm	st	st	pct	cs	ts	ctx"
print(header, file=out)
print(args[0],args[1])

def makeDict(fi1):
	dic1 = {}
	x,y = 0,1
	for i in fi1.readlines()[1:]:
		i = i.strip().split('\t')
		key1 = i[0]+":"+i[1]
		val1 = i[2:]
		val1[2] = int(val1[2])
		val1[3] = int(val1[3])
		dic1[key1] = val1
		x += 1
		if x >= 10000000:
			itr = y*10000000
			rtime = time.time()
			print("Make Dic %s lines Done, takes %s secs" % (itr,int(rtime-ctime)))
			x = 0 
			y += 1
	rtime = time.time()
	print("make dic Done, %s secs" % int(rtime-ctime))
	return dic1


dic1 = makeDict(fi1)
dic2 = makeDict(fi2)

x, y = 0, 1
for i,j in dic1.items():
	try: 
		dic2[i][2] += j[2]
		dic2[i][3] += j[3]
	except:
		dic2[i] = j
	x += 1
	if x >= 1000000:
		itr = y * 1000000
		rtime = time.time()
		print("merging %s lines done, takes %s secs" % (itr, int(rtime-ctime)))
		x = 0
		y += 1
rtime = time.time()
dic1 = []
key = sorted(dic2.keys())
print("merge done, #C is %s, takes %s secs" % (len(key),int(rtime-ctime)))

a, b = 0, 1
for k in key:
	x = k.split(':')
	print(x[0],x[1],sep="\t",end="\t",file=out)
	y = dic2[k]
	print(y[0],y[1],y[2],y[3],y[4],sep='\t',file=out)
	a += 1
	if a >= 10000000:
		itr = b * 10000000
		rtime = time.time()
		print("print %s lines done, takes %s secs" % (itr, int(rtime-ctime)))
		a = 0
		b += 1
out.close()
dic2 = []
rtime = time.time()
print("Print file done, takes %s secs" % int(rtime-ctime))
