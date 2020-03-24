#!/usr/bin/env python3
"""
Read FASTQ file & count total seq
multiple files can be applied and sum of them will output 

USAGE 
python3 f3_fq-count.py [fq] [fq.gz] [fq.bz2]

--LOG--
JAN.07th 2017	Change the file loading way, accept both gzip & bz2
NOV.07th 2015	Edit output & add average read length
OCT.23th 2015	Add python document
Oct.10th 2014	First build
"""


import sys, numpy, gzip
import bz2
"""
try:
	fi = sys.stdin
except:
	print(__doc__)
	sys.exit(2)
"""

def myopen(infile, mode="r"):
	if infile.endswith(".gz"):
		#print("GZIP file detected")
		return gzip.open(infile, mode=mode)
	elif infile.endswith(".bz2"):
		return bz2.open(infole, mode=mode)
	else:
		return open(infile, mode=mode)


fi = sys.argv
if len(fi) < 2:
	print(__doc__)
	sys.exit(2)
else:
	fi = fi[1:]
		

k,l = 0, 0
len_list = []
for f in fi:
	for n, m in enumerate(myopen(f)):
		if n % 4 == 1:
			m = m.strip()
			k += len(m)
			len_list.append(len(m))
			l += 1

print()
print("total nt    : %s" % k)
print("Num of reads: %s (if PE: %s)" % (l, int(l/2)))
print("Avg read len: %s" % round(k/l,2))
print("STDEV       : %s" % round(numpy.std(len_list),2))
print()
