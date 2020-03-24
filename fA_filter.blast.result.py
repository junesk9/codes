#!/usr/bin/env python3

"""
This tool is designed to filter the STACK resultant blast output from the RAD-seq
This tool perform to filter three steps;
	identity cutoff (default: 100)
	query/hit coverate cutoff (default: 100)
	pass only uniquely mapped result [when any query has two hits, the data is removed.]

The input file should be enhanced blast result format 6 from by f9_Enhance.Blast.Result.py 
& recommand to set value for "-max_target_seqs" option in blast >1

[USAGE]
fA_filter.blast.result.py [*enhanced.txt]

[output]
*.reduced.txt


JAN15th2016 Junesk9

"""


import sys


blast = sys.argv[-1]
prefix = ".".join(blast.split(".")[:-2])
blast = [i.strip().split("\t") for i in open(blast)]
header = blast[0]
blast = blast[1:]

perc_cutoff = 100
cov_cutoff = 100

if len(header) != 13:
	print("""
[WARNING] This code works well with "enhanced" blast result by f9_Enhance.Blast.Result.py 
[WARNING] While the input file seems not the proper data: %s
""" % prefix)
	sys.exit(2)
else: pass

print("[PROGRESS] Tottaly %s line submitted" % len(blast))
out = prefix + ".reduced.txt"
print("[PROGRESS] The output file %s is generated" % out)
out = open(out,"w")
print(*header,sep="\t",file=out)
filt = []
for i in blast:
	perc = i[2]
	cov = i[12]
	if float(perc) >= perc_cutoff and float(cov) >= cov_cutoff:
		filt.append(i)
	else: pass

name = [i[0] for i in filt]
foo = set()
common = [x for x in name if x in foo or foo.add(x)]
n_uniq = len(name)-len(common)
print("[PROGRESS] Filter %s lines; Uniq %s lines" % (len(name),n_uniq))
x,y = 0,1
for i in filt:
	if i[0] in common:
		pass
	else:
		print(*i,sep="\t",file=out)
		x += 1
	if x >= 10000:
		print("[PROGRESS] %s lines printed" % (y*10000))
		x = 0
		y += 1
print("[PROGRESS] %s lines printed, DONE" % n_uniq)
out.close()
	

