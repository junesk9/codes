#!/usr/bin/env python3
"""
Merge two 5mC input files to single file
Compare the 5mC site, add lines with different site to the 2nd input file
Equip sorting function

Usage:
python3 merge_inputs.py [1st input] [2nd input]

Feb. 2015 MoriToshi
"""
import sys

args = sys.argv[1:]
if not len(args) == 2:
	print("Only capable merge TWO files!")
	sys.exit()
prefix = args[0].split(".",1)[0]+".merged.input"
out = open(prefix,"w")


fi1 = open(args[0])
fi2 = open(args[1])

def ParseInput(fi):
	dic = {}
	for i in fi:
		i = i.strip()
		ch = i.split("\t")[0]
		site = i.split("\t")[1]
		key = ch+"-"+site
		dic[key] = i
	return dic
#end of def

dic1 = ParseInput(fi1)
dic2 = ParseInput(fi2)

key1 = set(dic1.keys())
key2 = set(dic2.keys())
diff = key1.difference(key2)
key1only = len(diff)
kommon = len(key1)-key1only
key2only = len(key2)-kommon
total1 = len(key1)+len(key2)-len(key1.intersection(key2)) 
total2 = key2only+kommon+key1only
print("total: %s, %s; file1: %s; file2: %s; common: %s" % (total1, total2, len(key1), len(key2), kommon))

print("#scf","site","str","mC","NmC","ctx","seq",sep="\t",file=out)
key3 = sorted(list(key2)+list(diff))
for i in key3:
	try:
		print(dic2[i], file=out)
	except:
		print(dic1[i], file=out)
#end of for
dic1,dic2,key3=[],[],[] ## Make computre memory empty
out.close()
