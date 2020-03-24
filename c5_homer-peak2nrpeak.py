#!/usr/bin/env python3

"""
Convert Homer FindPeak output "peaks.txt" to IGV-acceptable narrowPeak format.

generally combinate with for loop,
the Homer folder name given to the output file name

exp> for i in `find . -iname "peaks.txt"`; do tmp.homer-peak2nrpeak.py $i; done

Mar 1st 2017 Junesk9

"""


import sys

fi = sys.argv[-1]
prefix = fi.split("/")[-2]
fi = [i.strip() for i in open(fi)]
out = prefix+".homer-peak.narrowPeak"

out = open(out,"w")
for i in fi:
    if i.startswith("#"):
        if i.startswith("#Peak"):
            header=i.split("\t")
            header = ["#Chr"] + header[2:4] + ["Peak_ID","Normalized_tag_count"] + header[5:]
            print(*header,sep="\t",file=out)
        else:
            pass
    else:
        cnt = i.split("\t")
        cnt[0] = ":".join(cnt[0].split("-"))
        cnt = cnt[1:4] + [cnt[0],cnt[5]] + ["."] + cnt[6:]
        print(*cnt,sep="\t",file=out)

out.close()
            
