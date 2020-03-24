#!/usr/bin/env python3

"""
Collect only fasta listed in fasta_name.file
from a multi-fasta file
write on a new multi-fasta file


USAGE:

python3 fa_select.py -in [in_fa] -ls [in_ls] 

other options:
-size [intger]: only sequence length >= [size] will be retained [default: 0]
-except       : Only sequence not listed in the [in_ls] will be retaind [default: off]
-count        : stop collecting when the count leach to the given number [default: 0 (unlimited)]

===============Change log
Jan.27th.2017   Add "-count" option. only given number of seqs will outputs from top of the given list
Jan.18th.2017   Minor change (the gene name now become capitalized)
Oct.23rd.2105   Re-written for better readability
Aug.31th.2015   Change the Arguments parts
                Add new options "-no" "-size"
Dec.2nd 2014    Equip rapid data load procedure

"""

import sys

##################### Argument

opt = {"-in":"", "-ls":"","-size":0,"-except":False, "-count":0}
args = sys.argv[1:]

try:
    idx_in = args.index("-in")+1
    opt["-in"] = args[idx_in]
    idx_ls = args.index("-ls")+1
    opt["-ls"] = args[idx_ls]
except:
    print(__doc__)
    #print("USAGE: python3 f2_fa_select.py -in [] -ls []")
    sys.exit(1)

if "-except" in args:
    opt["-except"] = True
if "-out" in args:
    idx_out = args.index("-out")+1
    opt["-out"] = args[idx_out]
else:
    out = opt["-in"].split(".")[0]+"."+opt["-ls"].split(".")[0]+".out.fa"
    opt["-out"] = out
if "-count" in args:
    idx_cnt = args.index("-count")+1
    opt["-count"] = args[idx_cnt]

in_fa, in_ls, out, size, cnt = opt["-in"], opt["-ls"], opt["-out"],opt["-size"], int(opt["-count"])

print("""

    === Fasta multi-seq selector ===

    Running with
        Input fasta : %s
        Sequece list: %s
        Ouput file  : %s 
        Length limit: >=%s-nt [Default: 0]
        count mode  : %s [Default: 0 (unlimited)]
        Except mode : %s [Default: off]
""" % (in_fa, in_ls, out, size, cnt, opt["-except"]))

######################### Methods

def ParseFasta(in_fa):

    fa_dic = {}
    for i in open(in_fa):
        i = i.strip()
        if i.startswith(">"):
            name = i.split()[0].upper()[1:]
            fa_dic[name] = []
        else:
            fa_dic[name].append(i)
    for key in fa_dic.keys():
        fa_dic[key] = "".join(fa_dic[key])
        if fa_dic[key].endswith("."):
            fa_dic[key] = fa_dic[key][:-1]

    return fa_dic


########################## Main

fa_dic = ParseFasta(in_fa)
in_ls = [i.strip().split()[0].upper() for i in open(in_ls)]
in_ls = list(set(in_ls))
fa_ls = list(fa_dic.keys())

print(in_ls[:5])
print(fa_ls[:5])
#a = True if in_ls[0] in fa_ls else False
#print(a)


#sys.exit()

"""
fa_tmp = []
for i in fa_ls:
    i = i.split(".")[0]
    fa_tmp.append(i)
fa_ls = fa_tmp
print(in_ls[:10])
print(fa_ls[:10])
sys.exit()
"""

before = len(fa_ls)
err_count = 0

if opt["-except"] == True:
    for i in in_ls:
        try:
            fa_ls.remove(i)
        except ValueError:
            i = ">"+i
            try:
                fa_ls.remove(i)
            except ValueError:
                i = i.lower()
                try:
                    fa_ls.remove(i)
                except ValueError:
                    i = i[1:]
                    try:
                        fa_ls.remove(i)
                    except ValueError:             
                        print("Cannot find %s in the INPUT file, re-check the data" % i)
                        print(fa_ls)
                        sys.exit(1)
    after = len(fa_ls)
    print("[Except mode turned on] From %s seqs, remove %s seqs, remains %s seqs" % (before, len(in_ls) ,after))
elif opt["-except"] == False:
    tmp_ls = []
    for i in in_ls:
        if i.startswith(""):
            i = i[1:]
        if not "." in i:
            i = i + ".1" 
        try:
            seq = fa_dic[i]
            tmp_ls.append(i)
        except KeyError:
            err_count += 1
            print("[WARNING] Cannot find %s %s/%s " % (i, err_count, len(in_ls)))
    fa_ls = tmp_ls
    after = len(fa_ls)
    print("[Except mode turned off] From %s seqs, select %s seqs, obtains %s seqs"% (before, len(in_ls), after))

if not size == 0:
    cut_count = 0
    print("[Size cut mode] The seqs len <%s nt will be removed" % size)
    for i in fa_ls:
        if len(fa_dic[i]) < size:
            fa_ls.remove(i)
            cut_count += 1
    print("[Size cut mode] cut %s seq, remains %s seqs" % (cut_count, len(fa_ls)))

out = open(out,"w")
if cnt == 0 or cnt > len(fa_ls):
    for i in fa_ls:
        tag = ">" + i
        print(tag,file=out)
        print(fa_dic[i],file=out)
elif cnt < len(fa_ls):
    print("[Counted-output mode] First %s seqs will be printed" % cnt)
    for en,i in enumerate(fa_ls):
        if cnt > en:
            print(i,file=out)
            print(fa_dic[i],file=out)
        else:
            pass
out.close()
          


