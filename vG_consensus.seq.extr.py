#!/usr/bin/env python3

"""
Extract consensus seq from reference genome by parsing both GFF & VCF
 

[USAGE]
python3 vG_consensus.seq.extr.py -fa [genome_fasta] -gff [GFF3] -vcf [GATK-vcf] -id [gene_id] -mode [gene|cds|promoter] -out [out.file]

[output]
multifasta files each labeled >gene_id + accs

[options]
**mendatory options**
-fa  : multifasta of the reference genome
-gff : GFF3 file corresponing to "fa"
-vcf : merged gVCF resultant of GATK & requires ID column, [ch:pos:ref:alt] which can be added by following command
       [bcftools norm -Ou -m -any *.vcf | bcftools norm -Ob -f ref.fa | bcftools annotate -Ov -x ID -I +'%CHROM:%POS:%REF:%ALT' > *output]
-id  : gene id or specific region [ch:st-ed]
**other options**
-mode : [GENE|CDS|PROMOTER] GENE mode extracts whole mRNA region, or not genic region if specific region given from "-id"
                            CDS mode discards all introns and UTRs, only activated when gene id specified "-id"
-out  : specify the output file name
-indel: only for CDS mode. turn on the display intron/indels (appeared as */., respectibely) by claim this option.
-snpq : Choose vcf-SNP line by the variant_quality (4th column) to apply (default: 30)

May 2016 Junesk9

==== Change log ====
Jan 31th 2020	 Add SNP_quality option (-snpq) to chose SNP to apply
Dec 12th 2016	 Add optional view of VCF lines processed (default: ON)
Aug 10th 2016	 Add optional annotation steps for VCF files has not been annotated.
Jul 25th 2016	 now the strand of genome region input is evaluated by number of st:ed,
                 if st number is bigger than ed, the reverse strand will output.
Jun 14th 2016	 Fix bugs; add comment to pydoc
May 29th 2016    First  build
"""

import sys, time
try:
    import pandas as pd #ver0.18.1
except:
    print("""
    requires Pandas package (>=v0.18.1)
    try "pip3 install pandas --upgrade"
""")
    sys.exit()


################ survey Argument

args = {"fa":False ,"gff":False, "vcf":False, "mode":"GENE", "p_len":2000, "id ":False, "snpq":30, "disp_indel": False, "prefix": False, "vcfview":True}
## Available modes, "gene", "CDS", "promoter"
## -gff ** -vcf ** -mode ** [gene_id or range (chr:st-ed)]

"""
args = {"fa":"Arabidopsis_thaliana.TAIR10.31.dna.chromosome.3.fa",
        "gff":"tair10_m.gff", 
        "vcf":"AT.E180plus.chr3.151-152.reEval.vcf.recode.vcf",
        "mode":"GENE",
        "p_len":2000,
        "id":"3:15178872-15179100",
        "disp_indel": True,
        "prefix": False}
"""

arg = sys.argv[1:]
if len(arg) == 0:
    print(__doc__)
    sys.exit()
if "-gff" in arg:
    gff_idx = arg.index("-gff")+1
    args["gff"] = arg[gff_idx]
else: pass
#end of if
if "-vcf" in arg:
    vcf_idx = arg.index("-vcf")+1
    args["vcf"] = arg[vcf_idx]
else: pass
#end of if
if "-out" in arg:
    out_idx = arg.index("-out")+1
    args["prefix"] = arg[out_idx]
else: pass
#end of if
if "-fa" in arg:
    fa_idx = arg.index("-fa")+1
    args["fa"] = arg[fa_idx]
else: pass
#end of if
if "-indel" in arg:
    args["disp_indel"] = True
else: pass
#end of if
if "-id" in arg:
    id_idx = arg.index("-id")+1
    args["id"] = arg[id_idx].upper()
else: pass
#end of if
if "-mode" in arg: 
    mode_idx = arg.index("-mode")+1
    args["mode"] = arg[mode_idx].upper()
    if args["mode"] != "PROMOTER":
        args["p_len"] = False
    else: pass
else: pass
#end of if  
if "-snpq" in arg:
    snpq_idx = arg.index("-snpq")+1
    args["snpq"] = arg[snpq_idx]
else: pass
#end of if 



############### Display running options
fa,gff,vcf,mode,p_len,g_id,disp_indel = args["fa"],args["gff"],args["vcf"],args["mode"],args["p_len"],args["id"],args["disp_indel"]
prefix = args["prefix"]
if prefix == False:
    time = str(int(time.time()))[6:]
    if disp_indel == True and mode == "CDS":
        prefix = g_id.split(".")[0]+"."+mode+"."+time+".indel.consensus.fa"
    else:
        prefix = g_id.split(".")[0]+"."+mode+"."+time+".consensus.fa"
else: pass


print("""

        ### Consensus Sequence Extractor ###
        ###   via refereing GFF3 & VCF   ###

        GENE_ID    : %s (AGI code, or ch:st-ed, mendatory)
        run_mode   : %s (GENE, CDS, PROMOTER default: GENE)
        indel_disp : %s (display indel"."/intron"*", default: False)
    
        FA         : %s (mendatory)
        VCF        : %s (mendatory)
        GFF        : %s (default: False)
        promoter   : %s bps (default: 2,000)

        out_file   : %s

""" % (g_id, mode, disp_indel, fa, vcf, gff, p_len, prefix))


##########option validation 
if args["fa"] == False or args["vcf"] == False or args["id"] == False:
    print("""
[ERROR] all -vcf, -fa & -id options are mandetory!""")
    print(__doc__)
    sys.exit(1)
else: pass
#end of if
if not args["mode"] in ["GENE","PROMOTER","CDS"]:
    print("""
    [ERROR] mistypes in mode selecting, you typed: %s
    """ % args["mode"])
    sys.exit(1)

##################### Fasta parsing
fa_dic = {}
for i in open(fa):
    i = i.strip().split()[0].split("|")[0] ##discard unexpected things
    if i.startswith(">"):
        name = i[1:].upper() ## retrive chromosome identifier without ">"
        fa_dic[name] = []
    else:
        fa_dic[name].append(i)
    #end of if
for name in fa_dic.keys():
    fa_dic[name] = "".join(fa_dic[name])
#end of for
fa = fa_dic

### parsed fasta retrieved "fa_dic" 
print("[PROGRESS] Load Fasta done")
################### GFF parsing
gene_ls = {}
cds_ls = {}
cds_pos = {}
if gff != False:
    gff = [i.strip() for i in open(gff)]
    for line in gff:
        if not line.startswith("#"):
            line = line.split("\t")
            if line[2] in ["mRNA","transcript"]:
                gene = line[8].split(";")[0].split("=")[1] #aim to retrive AGI code from "ID="
                ch = line[0]
                st = line[3]
                ed = line[4]
                strn = line[6]
                gene_ls[gene] = [ch,st,ed,strn]
            elif line[2] == "CDS" and mode == "CDS":
                cds = line[8].split(";")[0].split("=")[1] #aim to retrive AGI code from "Parent="
                strn = line[6]
                st = int(line[3])
                ed = int(line[4])
                ch = line[0]
                exon_rng = list(range(st,ed+1))
                try:
                    cds_pos[cds] += exon_rng
                except KeyError:
                    cds_pos[cds] = exon_rng
                #end of try
                if strn == "+":
                    cds_min = min(cds_pos[cds])
                    cds_max = max(cds_pos[cds]) +2 ## I still have no idea, though when fwd-strand, only 2-bp is spared for terminater codon
                elif strn == "-":
                    cds_min = min(cds_pos[cds]) -3
                    cds_max = max(cds_pos[cds])
                else: pass
                #cds_min, cds_max = min(cds_pos[cds]), max(cds_pos[cds]) 
                #cds_min = min(cds_pos[cds]) - 3 if strn == "-" else min(cds_pos[cds])
                #cds_max = max(cds_pos[cds]) - 1 + 3 if strn == "+" else max(cds_pos[cds])
                cds_ls[cds] = [ch, cds_min, cds_max, strn]
            else: pass
            #end of if
        else: pass
        #end of if
    #end of for
else: pass
#end of if
# gene_ls, cds_ls, cds_pos are retrived
print("[PROGRESS] Load GFF done")
################### set range to survey

#rgn = str(arg[-1]).upper()
rgn = g_id
if ":" in rgn and "-" in rgn: #surveying region can be set as chr:st-ed
    gff = False ## The gff option would be ignored
    mode = "GENE" ## Prevent the misinput, basically "mode" is ignored with the genic region input.
    rgn_ch = rgn.split(":")[0]
    st = int(rgn.split(":")[1].split("-")[0])
    ed = int(rgn.split(":")[1].split("-")[1])
    if st > ed:
        rgn_st, rgn_ed = ed-1, st
        rgn_str = "-"
    else:
        rgn_st, rgn_ed = st-1, ed
        rgn_str = "+"
    """
    rgn_st = int(rgn.split(":")[1].split("-")[0])-1
    rgn_ed = int(rgn.split(":")[1].split("-")[1])
    rgn_str = "+"
    """
elif rgn.startswith("AT"): ## Gene id is expected, temporary, only AGI is applicapable now.
    rgn = rgn+".1" if not "." in rgn else rgn
    #print(rgn)
    if mode in ["GENE","PROMOTER"]:
        rgn = gene_ls[rgn]
    elif mode == "CDS":
        cds_pos = cds_pos[rgn]
        rgn = cds_ls[rgn]
    else: pass
    #end of if
    rgn_ch = rgn[0]
    rgn_st = int(rgn[1])-1
    rgn_ed = int(rgn[2])
    rgn_str = rgn[3]
    if mode == "PROMOTER":
        p_len = int(p_len)
        if rgn_str == "+":
            rgn_ed = rgn_st
            rgn_st = rgn_ed - p_len
        elif rgn_str == "-":
            rgn_st = rgn_ed
            rgn_ed += p_len
        else: pass
        #end of if
    else: pass
    #end of if
else: pass
#end of if
print(mode, rgn_ch, rgn_st, rgn_ed, rgn_str, sep="\t")
print(sorted(cds_pos)[:5])
print(sorted(cds_pos)[-5:])
#sys.exit(1)

###################### parsing a fasta based on rgn
fa = fa_dic[rgn_ch][rgn_st:rgn_ed]
#print(fa)
#sys.exit()
fa = [i for i in fa]

if mode == "CDS":
    cds_rng = list(range(rgn_st,rgn_ed+1))
    intron_pos = list(set(cds_rng) - set(cds_pos))
    #print(rgn_st, rgn_ed)
    #print(len(cds_rng),len(cds_pos),len(intron_pos))
    #sys.exit()
    for pos in intron_pos:
        idx = pos - rgn_st -1
        try:
            if not idx in [0,1,2,len(fa)-1,len(fa)-2,len(fa)-3,-1]: # CDS do not cover terminater codon
                fa[idx] = "*"
        except IndexError:
            print(len(fa), idx)
else: pass
#end of if

fa = "".join(fa)
#check
#print(fa)
#sys.exit()
#print(len(fa),fa[:10],fa[-10:])

            
###################### Parsing vcf based on rgn and dataframing using Pandas
tmp = []
for line in open(vcf):
    if not line.startswith("##"):
        line = line.strip().split()
        if line[0].startswith("#CHR"):
            if args["vcfview"] == True:
                print(*line,sep="\t")
            header = line
        elif line[6] == "PASS" or float(line[5]) >= args["snpq"]: ## collect the vcf data associated in RGN
            ch = line[0]
            site = int(line[1])
            if ch == rgn_ch and rgn_st-5 <= site <= rgn_ed+5:
                if args["vcfview"] == True:
                    print(*line,sep="\t")
                tmp.append(line)
            else: pass
        else: pass
        #end of if
    else: pass
    #end of if
#end of for
    
vcf = pd.DataFrame(tmp, columns=header)
#print(vcf["POS"].head(10))
vcf = vcf.sort_values(by="POS",ascending=True) # sorting by position
#print(vcf["POS"].head(10))
tmp = [] # blanks the memory
## print(vcf.shape)
## Retrieve vcf & header
if vcf.shape[0] < 1:
    print("[ERROR] no polymorphysm matched from VCF!")
    sys.exit(2)
else:
    print("[PROGRESS] Dataframing VCF done, %s polymorphysms retrieved" % vcf.shape[0])
#end of if

############ Functions
def revcomp(fa):
    tmp = ""
    for i in fa[::-1]:
        if i.upper() == "A":
            tmp += "T"
        elif i.upper() == "T":
            tmp += "A"
        elif i.upper() == "C":
            tmp += "G"
        elif i.upper() == "G":
            tmp += "C"
        elif i.upper() == "N":
            tmp += "N"
        elif disp_indel == True:
            tmp += i
        else: pass

    return tmp


###################### Extract consensus seq by accessions
header = header[9:]
print("[PROGRESS] Start extract consensus seq of %s accessions" % len(header))

vcf_pos = vcf["POS"].tolist()
#print(vcf_site[:5],vcf_site[-5:])


out = open(prefix,"w")
print("[PROGRESS] extracting TAIR10 reference ...")
print(">TAIR10"+"_"+g_id,file=out)
if disp_indel == False:
    ref_fa = "".join(i for i in fa if i not in ".*")
elif disp_indel == True:
    ref_fa = "".join(i for i in fa)
ref_fa = revcomp(ref_fa) if rgn_str == "-" else ref_fa
print(ref_fa, file=out)

for head in header:
    print("[PROGRESS] extracting %s ..." % head)
    tmp_indel = []
    tmp_fa = [i for i in fa]
    #print("".join(tmp_fa))
    #sys.exit()
    for site in range(rgn_st,rgn_ed+1):
        site = str(site)
        if site in vcf_pos:
            val = vcf.query("POS == @site")[head].tolist()[-1]
            #print(snp)
            #sys.exit()
            val = val.split(":")[0]
            if val == "0/0" or val == "./." or val =="0|0" or val == ".|.":
                pass
            else:
                val_id = vcf.query("POS == @site")["ID"].tolist()[-1]
                if val_id == ".": ## When the VCF has not been annotated AUG.10th.2016
                    val_site = vcf.query("POS == @site")["POS"].tolist()[-1]
                    val_ref = vcf.query("POS == @site")["REF"].tolist()[-1]
                    val_alt = vcf.query("POS == @site")["ALT"].tolist()[-1]
                else:
                    val_id = val_id.split(":")
                    val_site = val_id[1]
                    val_ref = val_id[2]
                    val_alt = val_id[3]
                #"""
                while len(val_ref) > len(val_alt): ## mark deletion ATG:A > ATG:A..
                    val_alt += "." 
                #end of while
                #### Actual substitution step
                #"""

                pos = int(val_site) - int(rgn_st) -1
                #print(val_ref)
                #print(val_alt)
                #print(fa[pos-1:pos+2])
                if fa[pos] not in ".*" and pos != -1: #and len(val_ref) == len(val_alt) == 1:
                    if fa[pos] != val_ref[0]:
                        print(fa)
                        print(val_id,val_site, rgn_st)
                        print("error fa: %s val: %s>%s unmatched" % (fa[pos-1:pos+2],val_ref,val_alt))
                        sys.exit()
                    val_alt = [i for i in val_alt]
                    try: 
                        tmp_fa[pos:pos+len(val_ref)] = val_alt
                    except:
                        print(pos,len(val_ref))
                        tmp_fa[3064]
                        print(val_alt)
                        sys.exit()
                else: pass
                #end of if
    
    if disp_indel == False:
        proc_fa = "".join(i for i in tmp_fa if i not in ".*")
    elif disp_indel == True:
        proc_fa = "".join(i for i in tmp_fa)
    #end of if
    #print(proc_fa)
    #sys.exit()
    proc_fa = revcomp(proc_fa) if rgn_str == "-" else proc_fa       
    print(">"+head+"_"+g_id,file=out)
    print(proc_fa,file=out)
out.close()
print("[PROGRESS] All progress has done")
sys.exit(2)

