#!/usr/bin/env python3
"""
ADD REFERENCE data to PLINK PED file by refering MAP file,

[USAGE]

./PED.add-ref.py [prefix].map >> [prefix].ped


FEB2016 Junesk9
"""


import sys

fi = sys.argv[-1]

db = ["H06_Col-0","H06_Col-0",0,0,0,0]

for i in open(fi):
	i = i.split()[1]
	ref = i.split(":")[2]
	db.append(ref)
	db.append(ref)

print(*db,sep="\t")

