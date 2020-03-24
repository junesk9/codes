#!/usr/bin/env python3
"""
De-duplicator of multi-fasta file

compare both reverse complementary strand,
only unique sequences are retrived

This performs very fastly using a Counter packages


Jun 1st 2016 Junesk9

"""


import sys, time
#from hashlib import md5
from collections import Counter

fi = sys.argv[-1]


def revcomp(fa):
        tmp = ""
        for f in fa[::-1]:
                if f.upper() == "A":
                        tmp += "T"
                elif f.upper() == "T":
                        tmp += "A"
                elif f.upper() == "G":
                        tmp += "C"
                elif f.upper() == "C":
                        tmp += "G"
                else:
                        tmp += f
        return tmp

dedup = True
ctime = time.time()

if dedup == True:
        print("[PROGRESS] Start de-duplication")

        dic, dic_r = {},{}
        dic2 = {}
        x, n = 0, 0
        for i in open(fi):
                i = i.strip()
                if i.startswith(">"):
                        name = i
                        dic[name] = ""
                        dic_r[name] = ""
                        x += 1
                else:
                        dic[name] = i.upper()[:-1]
                        dic2[dic[name]] = name
                        dic_r[name] = revcomp(i.upper()[:-1])
                        #print(dic_rmd[name])
                        #sys.exit()
                if x > 99999:
                        x = 0
                        n += 100000
                        print("[PROGRESS] %s seqs are loaded, takes %s sec" % (n, int(time.time()-ctime)))	
        print("[PROGRESS] loading fasta done, %s seqs are loaded, takes %s min" % (n+x, int((time.time()-ctime)/60)))


        if len(dic.keys())*2 == len(set(list(dic.values())+list(dic_r.values()))):
                print("[DONE] No duplicate is found; All progress done")
        else:
                print("[PROGRESS] Duplication was detected, de-duplication step activated")
                #keys = list(dic.keys())
                #before = len(keys)
                #dup = []
                seqs = list(dic.values()) + list(dic_r.values())
                #print("[PROGRESS] start counting")
                count = Counter(seqs)
                #print("[PROGRESS] end counting, takes %s min" % int((time.time()-ctime)/60))
                y, z, m = 0,0,0
                prefix2 = ".".join(fi.split(".")[:-1])+".dedup.fa"
                out2 = open(prefix2, "w")
                for i,j in count.items():
                        y += 1
                        if j == 1:
                               try:
                                      name = dic2[i]
                                      print(name, i,sep="\n", file=out2)
                                      m += 1
                               except KeyError: pass
                        """
                        else: 
                                print(i,j,sep="\n")
                                sys.exit()
                        """
                        if y > 99999:
                                y = 0
                                z += 100000
                                print("[PROGRESS] %s seqs are processed, takes %s min" % (z, int((time.time()-ctime)/60)))
                out2.close()
                seqs, dic, dic_r, dic2 = [],[],[],[]
                print("[DONE] Unique %s seqs retained from %s seqs" % (m, n+x))
                """
                x, y, z = 0,0,0
                for key in keys:
                        if seqs.count(dic[key]) > 1:
                                dup.append(key)
                        else: pass
                        y += 1
                        if y > 499:
                                y = 0
                                z += 500
                                print("[PROGRESS] %s seqs are processed, takes %s min" % (z, int((time.time()-ctime)/60))
                keys = keys - dup
                after = len(keys)
                print("[PROGRESS] Deduplication done, %s to %s" % (before,after))
                print("[PROGRESS] Start writing the output file")
                prefix2 = ".".join(prefix.split(".")[:-1])+".dedup2.fa"
                out2 = open(prefix2, "w")
                for key in keys:
                        print(key, dic[key],sep="\n", file=out2)
                out2.close()
                """
                #print("[DONE] the output file % genereated, all process has done" % prefix2)
