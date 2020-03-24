#!/usr/bin/env python3

"""
Survey Chromosome distribution of assigned DMCs &DMRs
input DMR file is as output from DMC.DMR analysis
Chromosone length & other parameter is set for A. thaliana now
(Change [ch_len] [win] for other organisms and different window size)

Usage:
python3 m8_dmr.chr.dist.py [*DMC] # for DMC
python3 m8_dmr.chr.dist.py [*hyper.DMC] [*hypo.DMC] # for DMR

Oct 21st 2014 by ToshiMori

---- Log ----
Nov.11.14 Debugging for usage of DMC

Nov.01.14 Automatically detect whether input is DMC or DMR by input file number (1 or 2)
           Also automatically parse two DMR files properly by its file name contains (hypo/hyper)
           And Now make single output file contains both hyper/hypo count info for each seq context
           (Six outputs -> Three outputs, more convinient for adapt to R)

Oct.25.14 Also accept DMC file as input; verify DMR or DMC input autometically (by counting column no.)
          Alternativly adapt (-) to count, when the input file contain 'hypo'.

"""


import sys

########################################################
########### Args

ch_len = {1:30427671, 2:19698289, 3:23459830, 4:18585056, 5:26975502} # In case of At
win = 500000 # size of single window (bp)
fi = sys.argv[1:]
prefix = ".".join(fi[0].split(".")[:2])

## Check input form
if len(fi) > 2 or len(fi) == 0:
    print("""
    ##### Input file error ####
    The No. input file shoud be 1 for DMC
    or 2 for DMR (hyper/hypo)

    Nov.1st 2014 Toshimori
        """)

########################################################
########## Method

def PrepList():
    CG = []
    for i in sorted(ch_len.keys()):
        num = int(ch_len[i]/win)+1
        tmp = [0 for i in range(num)]
        CG.append(tmp)
    return CG

def PrintResult(CG, CGm, ctx):
    if len(sys.argv[1:]) == 1:
        name = prefix+"_"+ctx+"_DMC.chr.dist.txt"
    elif len(sys.argv[1:]) == 2:
        name = prefix+"_"+ctx+"_DMR.chr.dist.txt"
    out = open(name, 'w')

    # Count maximum length of chrs
    tmp = []
    for i in CG:
        tmp.append(len(i))
    max_num = max(tmp)
    
    print("ch1","ch2","ch3","ch4","ch5","ch1(-)","ch2(-)","ch3(-)","ch4(-)","ch5(-)",sep='\t',file=out)
    for i in range(max_num):
        try: print(CG[0][i], end="\t", file=out)
        except: print('NA', end="\t", file=out)
        try: print(CG[1][i], end="\t", file=out)
        except: print('NA', end="\t", file=out)
        try: print(CG[2][i], end="\t", file=out)
        except: print('NA', end="\t", file=out)
        try: print(CG[3][i], end="\t", file=out)
        except: print('NA', end="\t", file=out)
        try: print(CG[4][i], end="\t", file=out)
        except: print('NA', end="\t", file=out)
        try: print(CGm[0][i], end="\t", file=out)
        except: print('NA', end="\t", file=out)
        try: print(CGm[1][i], end="\t", file=out)
        except: print('NA', end="\t", file=out)
        try: print(CGm[2][i], end="\t", file=out)
        except: print('NA', end="\t", file=out)
        try: print(CGm[3][i], end="\t", file=out)
        except: print('NA', end="\t", file=out)
        try: print(CGm[4][i], end="\t", file=out)
        except: print('NA', end="\t", file=out)
        print(file=out)
    out.close()

######################################################
########## Main

ch = sorted(ch_len.keys())
CG = PrepList()
CHG = PrepList()
CHH = PrepList()
CGm = PrepList()
CHGm = PrepList()
CHHm = PrepList()

if len(fi) == 1: # A DMC file contains both hyper/hypo info.
    fi1 = fi[0]
    for i in open(fi1).readlines()[1:]:
        i = i.strip().split('\t')
        ctx = i[-1].split(',')[0]
        tmp = [ctx, int(i[0]), int(i[1]), float(i[6])]
        loc = int(tmp[2]/win)
        if ctx == "CG" and tmp[3] < 0:
            CGm[ch.index(tmp[1])][loc] -= 1
        elif ctx == "CHG" and tmp[3] < 0:
            CHGm[ch.index(tmp[1])][loc] -= 1
        elif ctx == "CHH" and tmp[3] < 0:
            CHHm[ch.index(tmp[1])][loc] -= 1
        elif ctx == "CG" and tmp[3] > 0:
            CG[ch.index(tmp[1])][loc] += 1
        elif ctx == "CHG" and tmp[3] > 0:
            CHG[ch.index(tmp[1])][loc] += 1
        elif ctx == "CHH" and tmp[3] > 0:
            CHH[ch.index(tmp[1])][loc] += 1
elif len(fi) == 2: # DMR input have 2 files (hyper/hypo)
    fi1 = sorted(fi) # Make file order as hyper -> hypo
    for n, m in enumerate(fi1):
        for i in open(m).readlines()[1:]:
            i = i.strip().split('\t')
            ctx = i[-1].split(',')
            for j in ctx:
                tmp = [j, int(i[1]), int(i[2])]
                loc = int(tmp[2]/win)
                if j == "CG" and n == 1:
                    CGm[ch.index(tmp[1])][loc] -= 1
                elif j == "CHG" and n == 1:
                    CHGm[ch.index(tmp[1])][loc] -= 1
                elif j == "CHH" and n == 1:
                    CHHm[ch.index(tmp[1])][loc] -= 1
                elif j == "CG" and m == 0:
                    CG[ch.index(tmp[1])][loc] += 1
                elif j == "CHG" and m == 0:
                    CHG[ch.index(tmp[1])][loc] += 1
                elif j == "CHH" and m == 0:
                    CHH[ch.index(tmp[1])][loc] += 1        
  
PrintResult(CG, CGm, "CG")
PrintResult(CHG, CHGm, "CHG")
PrintResult(CHH, CHHm, "CHH")

