#!/usr/bin/env python3
"""
Find proximal genes around a DMR
the proximal distanct = win


Usage:
python3 find.DMR.prox.gene.py [gff] [DMR]
(The DMR file should contain from 1st column; [DMR_name][chr][st][ed])

The output appears as STDOUT

Oct. 2014 by ToshiMori

=========change log
15.02.12 Also adaptable with "s5_bedPeakMaker.py" output
14.10.28 Also print context at the end of each line
"""
import sys

###############################################################
############## Arguments

args = sys.argv[1:]

fi = open(args[1])
gff = open(args[0])
#out = open('hyper.DMR.gene.txt','w')
win = 2000 # proximal distances from the borders of a DMR

################################################################
############# Body

dic_fi = {}
for i in fi.readlines()[1:]:
    i = i.rstrip().split('\t')
    dic_fi[i[0]] = i[1:]
ln = len(dic_fi)

dic_gff = {}
dic_gff2 = {}
for i in gff.readlines():
    i = i.strip().split('\t')
    name = i[-1].split(';')[0].split('=')[-1]
    dic_gff[name] = [i[0], i[3], i[4]]
    dic_gff2[name] = [i[2], i[6]]

x,y = 0, 1
count = 0
print("#DMR", "status", "rel_dist", "Gene", "class", "chr","str","gene_st","gene_len","ctx(depth)" ,sep="\t") 
for j,k in dic_fi.items():
    for l,m in dic_gff.items():
        if k[0]+"|" in m[0] or m[0] == k[0]:
            if int(k[1])-win <= int(m[1]) <= int(k[2])+win or int(k[1])-win <= int(m[2]) <= int(k[2])+win:
                rng1 = {i for i in range(int(k[1])-win,int(k[2])+win+1)}
                rng2 = {i for i in range(int(m[1]),int(m[2])+1)}
                comm = sorted(list(rng1.intersection(rng2)))
                gene_len = int(m[2])-int(m[1])+1
                if len(comm) != 0:
                    if comm[-1] < int(k[1]):
                        status = "upstream"
                        dist = comm[-1] - int(k[1])
                    elif comm[0] > int(k[2]):
                        status = "downstream"
                        dist = comm[0] - int(k[2])
                    else:
                        status = "overlapping"
                        dist = 0
                    print(j,status,dist,l,dic_gff2[l][0],k[0],dic_gff2[l][1],m[1],gene_len,k[-1], sep='\t')
                    count += 1
    x += 1
    if x*10 >= ln:
        print("%d%% Done, takes %d secs" % (y*10, int(time.time()-ctime)))
        x = 0
        y += 1

#out.close()
'''
print()
print("No. DMR: %d" % ln)
print("No. gene: %d" % count)
'''

