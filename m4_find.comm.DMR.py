#!/usr/bin/env python3
"""
Find Common DMR regions from two DMR out files
(fi1, fi2)
Output Common DMR files with seq-info & coverage perc
(stdout)

Usage:
python3 find.DMR.prox.gene.py [DMR1] [DMR2]
(The DMR file should contain from 1st column; [DMR_name][chr][st][ed])

Oct. 2014 by Toshimori 


Log ====
14.10.28 - Also output common seq context info & Eject DMR with no common context
14.10.28 - Also outpur Avg 5mC level for each context
14.10.28 - convert some 5mC values appeared as -100 to 'NA', for convinience in use of R.  those values have no corresponding context.
"""

import sys

#############################################################
############# Arguments

args = sys.argv[1:]
fi1 = open(args[0])
fi2 = open(args[1])
#out = open('comm.DMR.hyper.txt','w')

#############################################################
############ Methods

def ReadFi(fi1):    
    dic = {}
    for i in fi1.readlines()[1:]:
        i = i.strip().split('\t')
        dic[i[0]] = i[1:]
    return dic

############################################################
########### Body

dic_fi1 = ReadFi(fi1)
dic_fi2 = ReadFi(fi2)
ln1 = len(dic_fi1)
ln2 = len(dic_fi2)

x,y = 0, 1
count = 0
print("#Comm","Chr","Comm_st","Comm_ed","DMR1","DMR1_st","DMR1_ed","DMR2","DMR2_st","DMR2_ed","DMR1/Comm","DMR2/Comm","s1_CG","s2_CG","s1_CHG","s2_CHG","s1_CHH","s2_CHH","ctx", sep="\t")
for j,k in dic_fi1.items():
    for l,m in dic_fi2.items():
        if k[0]+"|" in m[0] or m[0]+"|" in k[0] or m[0] == k[0]: ## To detect variable formats of scaffold names
            if int(k[1]) <= int(m[1]) <= int(k[2]) or int(k[1]) <= int(m[2]) <= int(k[2]):
                rng1 = {i for i in range(int(k[1]),int(k[2])+1)}
                rng2 = {i for i in range(int(m[1]),int(m[2])+1)}
                comm = sorted(list(rng1.intersection(rng2)))
                ctx = set(m[-1].split(',')).intersection(set(k[-1].split(',')))
                ctx = ','.join(sorted(ctx))
                for o in k[5:11]: # Some avg 5mC level appears as -100, since no context detected at the area. change -100 to NA for convinience in use of R
                    if str(o) == "-100":
                        o = "NA"
                    else: pass
                if len(comm) != 0 and len(ctx) >= 2:
                    count += 1
                    ID = "Comm"+str(count)
                    rate1 = round(len(comm)/len(rng1)*100,2)
                    rate2 = round(len(comm)/len(rng2)*100,2)
                    print(ID,k[0],comm[0],comm[-1],j,k[1],k[2],l,m[1],m[2],rate1,rate2,k[5],k[8],k[6],k[9],k[7],k[10],ctx, sep='\t')
'''              
    x += 1
    if x*10 >= ln1:
        print("%d%% Done, takes %d secs" % (y*10, int(time.time()-ctime)))
        x = 0
        y += 1
out.close()

    
print()
print("No. DMR1: %d" % ln1)
print("No. DMR2: %d" % ln2)
print("No. Comm: %d" % count)
'''
                
