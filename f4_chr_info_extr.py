#!/usr/bin/env python3

'''
Import multifast file with 'cat'
and output as
'chr	len
 chr	len
 chr	len
 ...	   '

to adapt genomeCoverageBed parameter -g
'''


import sys

fi = open(sys.argv[-1])

dict = {}
for i in fi:
	if i.startswith('>'):
		key = i.split()[0].strip()[1:]
		dict[key] = 0
	else:
		i = i.strip()
		value = len(i)
		dict[key] += value

for j, k in dict.items():
	print(j,k, sep='\t')
