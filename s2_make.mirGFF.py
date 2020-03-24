#!/usr/bin/env python3
"""
Parse mir_candidate.txt to .gff3 file
Oct. 2014 by ToshiMori
"""


fi = open('miR_candidate.txt')
out = open('RsMIR_v0.2.gff','w')

for j,i in enumerate(fi):
        i = i.strip().split('\t')
        try:
            annt = "ID="+i[0]+";Count="+i[5]+";Family="+i[14]+":"+i[15]+";Evidence="+i[16]
        except IndexError:
            annt = "ID="+i[0]+";Count="+i[5]+";Evidence="+i[16]
        print(i[1],"rsa-MIR_v0.2","precursor",i[3],i[4],".",i[2],".",annt,sep='\t', file=out)           
        annt = "ID="+i[0]+"-5p;Parent="+i[0]+";Seq="+i[6]
        print(i[1],"rsa-MIR_v0.2","mature-5p",i[7],i[8],".",i[2],".",annt,sep="\t", file=out)
        annt = "ID="+i[0]+"-3p;Parent="+i[0]+";Seq="+i[10]
        print(i[1],"rsa-MIR_v0.2","mature-3p",i[11],i[12],".",i[2],".",annt,sep="\t", file=out)
out.close()
