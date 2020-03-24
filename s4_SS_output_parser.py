#!/usr/bin/env python3
"""
Parse ShortStack output 'MIRNAs' folder
to single output txt file

Nov. 2014 ToshiMori
"""
import glob

input_ls = sorted(glob.glob("MIRNAs/*txt"))
out_fi = open("SS_miR_out.txt",'w')

print("#candidate","chr","orientation","pre_st","pre_ed","pre_seq","5p_st","5p_ed","5p_seq","3p_st","3p_ed","3p_seq",sep='\t', file=out_fi)
for i in input_ls:
    name = i.split('/')[1].split(".")[0]
    for n,m in enumerate(open(i)):
        if n == 15:
            m1 = m.strip().split(" ")
            ori = m1[-1]
            ch = m1[4].split(':')[0]
            sted = m1[4].split(':')[1]
            st = int(sted.split('-')[0])
            ed = int(sted.split('-')[1])
        if n == 17:
            seq = m.strip()
        if m.startswith("Seq"):
            m3 = m.strip().split(" ")
            m5p = m3[1]
            m3p = m3[-1]
            loc_m5p = seq.index(m5p)
            loc_m3p = seq.index(m3p)
    if ori == "+":
        m5p_st = st+loc_m5p
        m5p_ed = m5p_st+len(m5p)-1
        m3p_st = st+loc_m3p
        m3p_ed = m5p_st+len(m3p)-1
    elif ori == "-":
        m5p_ed = ed-loc_m5p
        m5p_st = m5p_ed - len(m5p)+1
        m3p_ed = ed-loc_m3p
        m3p_st = m3p_ed - len(m5p)+1
    print(name,ch,ori,st,ed,seq,m5p_st,m5p_ed,m5p,m5p_st,m5p_ed,m3p,sep='\t', file=out_fi)
out_fi.close()

