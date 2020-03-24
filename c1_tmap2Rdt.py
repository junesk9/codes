#!/usr/bin/env python3

"""
Convert cuffcumpare resultant *tmap file to Rdata file, which annotating
Cufflink's TCONS_ nomenclatures to the reference Gene codes.

[USAGE] 
Python3 tmap2Rdt.py [*tmap]

[Output]
*.Rdata.txt


[LOG]
Oct 4th 2015 First build


"""


import sys
fi = sys.argv[-1] 

prefix = ".".join(fi.split(".")[:2])+".Rdata.txt"
out = open(prefix,"w")
print("t_name","class_code","AGI_code",sep="\t",file=out)

fi = open(fi)
for i in fi.readlines()[1:]:
    i = i.rstrip().split("\t")
    print(i[4],i[2],i[1],sep="\t",file=out)
out.close()

