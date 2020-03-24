# ~/usr/bin/env python3
"""
For bwa-meth results, *.bed file,
survey the read depth for each Cytosyin globally,
output the [read_depth_# portion]

Oct.2014

=============Change log
Nov.19.2014 change input form [sys.stdin => sys.argv] 
"""


import sys

fi = sys.argv[1]

dict = {}
tot = 0
for i in open(fi).readlines()[1:]:
	i = i.strip().split('\t')
	k = int(i[4])+int(i[5])
	tot += 1
	try:
		dict[k] += 1
	except:
		dict[k] = 1
dict_ls = [int(i) for i in dict.keys()]
dict_ls = sorted(dict_ls)

for j in dict_ls:
	print(j, end='\t')
	print(dict[j], end='\t')
	perc = round(int(dict[j])*100/tot, 2)
	print(str(perc)+' %')
