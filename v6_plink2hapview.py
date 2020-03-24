#!/usr/bin/env python3

"""
PLINK inputs (*ped/*map) converter for Haploview usage

To apply plink inputs to run in Haploview (v1.9 present), need to process two things;
1. A C G T in the PED file sould be re-assigned as 1 2 3 4, respectively.
   (And also convert *[single-base deletion, assigned by PLINK[v >1.9]] to 0[missed])
2. The input files should be separated by chromosomes, since HapView is not enough memory-efficient than PLINK

This tool processes two plink input files and generates two Haploview-adaptable modified files..


[USAGE]
python3 v6_plink2hapview.py [prefix of PLINK inputs, e.g) sample.map/sample.ped -> sample] -m [specific SNP id] -l [integer]
	-m: specify an marker, then program will elutes markers surrounding the specified markers in 250k (default, adjustable by "-l")
	-l: specify the range to survey markers from the specified marker by "-m" option, only activated when "-m" option is assigned.


[output]
---when no "-m" assigned
prefix.hapview.ch1.map/ped; prefix.hapview.ch2.map/ped; ....

---when "-m" assigned
prefix.hapview.[markername].map/ped


The output come suitable to adapt haploview, such as
>java -jar Haploview.jar -memory 200000 -n -pedfile [] -map [] -skipcheck -blockouput GAB -dprime

and the output is used for the LD block visualization, using Haploview without "-n: run without any visualization" or LD-plus
>java -jar Haploview.jar -pedfile [] -map []



=== Log ===
MAY 23 2016	Correct bugs, change filename generation
SEP 29 2015	Add "-l" option
SEP 25 2015	Add "-m" option
SEP 23 2015	First build

"""

import sys

################################### Args

args = sys.argv[1:]
if len(args) == 0:
    print(__doc__)
    sys.exit(1)

plink = args[0]
pedf = plink+".ped"
mapf = plink+".map"
marker_mode = False

marker, survey_rng = False, False
if "-m" in args:
    marker_mode = True
    idx_marker = args.index("-m")
    marker = args[idx_marker+1]
    marker_info = marker.split(":")
    marker_ch = marker_info[0]
    marker_bp = int(marker_info[1])
    survey_rng = "100000"
    if "-l" in args:
        idx_len = args.index("-l")
        survey_rng = args[idx_len+1]
        try:
            survey_rng = int(survey_rng)
        except:
            if survey_rng.endswith("k") or survey_rng.endswith("K"):
                survey_rng = int(survey_rng[:-1]+"000")
            elif survey_rng.endswith("m") or survey_rng.endswith("M"):
                survey_rng = int(survey_rng[:-1]+"000000")
            else:
                print("input range %s is invaild, only k/m is acceptable" % survey_rng)
                sys.exit(1)

mark_name = str(marker)
if ":" in marker: ## generally marker name as chr:bp:ref:alt, ":" is not good for the file name, so change to "chr-bp"
    mark_name = str("-".join(marker.split(":")[:2]))


if marker_mode == True:
	print("""

           This run is acheived with the marker-specific mode
		
		map       : %s
		ped       : %s
		marker    : %s
                survey_rng: %s (default: 100k+100k)

        """ % (mapf,pedf,marker,survey_rng))

#####################################Load MAP file and prepare outputs

if marker_mode == True:
    mapf = [i.strip().split("\t") for i in open(mapf)]

    #Find chromosome number "ch" that the SNP included
    ch = ""
    if ch == "":
        for i in mapf:
            map_ch = i[1].split(":")[0]
            if map_ch == marker_ch:
                ch = i[0]
            else: pass
    else: pass

    # Assess maximum bp of particular chromosome
    bps = []
    for i in mapf:
        if i[0] == ch:
            bps.append(int(i[-1]))
    max_bp = max(bps)

    # Assign the range of SNP data to collect, around 250 kbs
    rng_min = marker_bp - survey_rng
    if rng_min < 0:
        rng_min = 0
    rng_max = marker_bp + survey_rng
    if rng_max > max_bp:
        rng_max = max_bp
    rng = [rng_min, rng_max]

    ## Parse MAP file & obtain the SNP range for PED parsing
    en_rng = []
    map_t, map_k = 0,0
    outf = plink+"."+mark_name+".map"
    out = open(outf,"w")
    for en, i in enumerate(mapf):
        map_t += 1
        if i[0] == ch and rng[0] <= int(i[-1]) <= rng[-1]:
            print(*i, sep="\t",file=out)
            en_rng.append(en)
            map_k += 1
    en_max = max(en_rng)
    en_min = min(en_rng)
    en_rng = [en_min, en_max] ## Final info for PED parsing
    en_rng = [en_min*2-2, en_max*2-1]
    print("%s file generated, %s from %s snps selected, %s:%s - %s range" % (outf,map_k,map_t,ch,rng_min,rng_max)) 
    out.close()

    
    
elif marker_mode == False:
    map_sep = {} ## chr:No. SNP
    mapf = open(mapf)
    for i in mapf:
        i = i.strip()
        ch = i.split("\t")[0]
        try:
            map_sep[ch].append(i)
        except KeyError:
            map_sep[ch] = [i]



    chrs = sorted(list(map_sep.keys()))
    prefix = ".".join(plink.split(".")[:-1])
    for i in chrs:
        out = prefix+".hapview.ch"+str(i)+".map"
        out = open(out,"w")
        print(*map_sep[i],sep="\n",file=out)
        out.close()

    for i,j in map_sep.items():
        map_sep[i] = len(j)
    chr_len = map_sep


############################ Parse Convert string SNP to numeric SNP in PLINK resultant PED file

ped = [i.strip().split("\t") for i in open(pedf)]

## Separate and store header (column1-6)
header = []
for i in ped:
    #print(i[:10])
    tmp = i[:6]
    tmp = "\t".join(tmp)
    header.append(tmp)
    del i[:6]
#print(ped[0][:10])

## Convert ACGT to 1234	
for en,i in enumerate(ped):
    ind = header[en].split("\t")[0]
    for en,j in enumerate(i):
        if j == "A":
                i[en] = "1"
        elif j == "C":
                i[en] = "2"
        elif j == "G":
                i[en] = "3"
        elif j == "T":
                i[en] = "4"
        elif j == "*":
                i[en] = "0"
        elif not j == "0":
            print("[Warning] unxpected genotype! ind: %s colume: %s, appeares: %s" % (ind,en,j))
            sys.exit()

#print(len(ped[0]),sum(chr_len.values()))


######## IF marker selected
if marker_mode == True:
    en_rng = [en_min*2-2, en_max*2]
    ped_sel = []
    for i in ped:
        tmp = i[en_rng[0]:en_rng[-1]]
        tmp = "\t".join(tmp)
        ped_sel.append(tmp)

    outf = plink+"."+mark_name+".ped"
    out = open(outf,"w")
    for en, i in enumerate(header):
        line = i+"\t"+ped_sel[en]
        print(line,file=out)
    print("%s file generated" % outf)
    out.close()
#end of if


######## If marker not selected, parse file by chromosomes
elif marker_mode == False:
    ped_sep = []
    for en,i in enumerate(chrs):
        ped_sep.append([])
        nsnp = int(chr_len[i])*2
        for j in ped:
            #if en == 9: # Validate the parsing properity
                #print(nsnp, len(j))
                #sys.exit()
            tmp = j[:nsnp]
            tmp = "\t".join(tmp)
            del j[:nsnp]
            ped_sep[en].append(tmp)

    for en,i in enumerate(ped_sep):
        out = prefix+".hapview.ch"+str(en+1)+".ped"
        #print(len(i)) # should be the sample size
        #sys.exit()
        out = open(out,"w")
        for en2,j in enumerate(i):
            line = header[en2]+"\t"+j
            print(line,file=out)
        out.close()
            
