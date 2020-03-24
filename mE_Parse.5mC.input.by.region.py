#!/usr/bin/env python3


"""
Regional 5mC data collector from *.5mC.input files.


[USAGE]
Parse.5mC.input.by.region.py -p [prefix] 4:13019150-13019500 (ch:st-ed)

[OUTPUT]
A tab-separated plain txt file can be directly applied to boxplot() in R


Oct. 26th 2017 Junesk9

######## Log
FEB05. 2018	New option "-p" to set PREFIX of the output
                TIME score added to unique output file neames in everytime.
Dec22. 2017	Add context "CTX" colume data; and re-new the parsing body, mainly depending on Pandas

"""

import sys
import glob
import pandas 
import numpy
import time

################################################ Arguments

ctime = str(time.time())[-4:]
if len(sys.argv) == 1:
  print("[ERROR] Analyzing region should be specified !!")
  #print("[ERROR] Try agian: ch:st-ed")
  print(__doc__)
  sys.exit()
else:
  args = sys.argv
  if "-p" in args:
    p_idx = args.index("-p")+1
    prefix = args[p_idx]
    prefix = prefix
  else:
    prefix = "5mC-out"

  rgn = sys.argv[-1]
  try:
    ch = rgn.split(":")[0]
    se = rgn.split(":")[1].split("-")
    se = [int(se[0]),int(se[1])]
    st = min(se)
    ed = max(se)
  except:
    print("[ERROR] region specified should be like below: 1:10000-20000")
    print("[ERROR] Try agian: ch:st-ed")
    sys.exit()
print("[PROGRESS] Region specified like below %s:%s-%s" % (ch,st,ed))

inp_ls = glob.glob("*input")
inp_ls = sorted(inp_ls)
#inp_ls =["wt-minus.5mC.more4.input","rdd-z1.5mC.more4.input"]
print("[PROGRESS] Totally %s 5mC.input files are detected" % len(inp_ls))


################ Generate empty dataframe to fill after

#print(st, ed, ch)
#sys.exit()

columns = [i.split(".5mC")[0] for i in inp_ls]
columns = ["CTX"]+columns
index = [i for i in range(st,ed+1)]

df_ = pandas.DataFrame(index=index, columns=columns)
df_ = df_.fillna("NA") ## Use "NA" rather than "NaN" as NULL in R
#print(df_)
#sys.exit()


################ Fill the dataframe

for inp in inp_ls:
  inp = inp.strip()
  head = inp.split(".5mC")[0]

  print("[PROGRESS] Start parsing %s" % inp)
  for i in open(inp):
    i = i.strip().split()
    i_ch = i[0]
    i_pos = int(i[1])
    if i_ch == ch:
      ## pull info
      if st <= i_pos <= ed:
        mC = int(i[3])
        NmC = int(i[4])
        try:
          rate = round(mC/(NmC+mC),4)
        except TypeError:
          print(inp, i)
          sys.exit()
        i_ctx = i[5]
        
        #### refer CTX, if matched, fill the particular cell
        ctx = df_.ix[i_pos]['CTX'] 
        if ctx in ["NA", i_ctx]:
          df_.at[i_pos, 'CTX'] = i_ctx
          df_.at[i_pos, head] = rate
        else:
          print("[ERROR] GIVEN context did not matched; %s:%s %s vs. %s" % (i_ch,i_pos,ctx, i_ctx))
          df_.at[i_pos, head] = rate
          #sys.exit()
  #print(df_)
  #sys.exit()



############# Delete lines with "NA" for CTX

df2 = df_[df_["CTX"] != "NA"]
#print(df2)
#sys.exit()

############ Print the  output

outf = prefix + "." + ctime + "." + "parse.txt"
#outf = "5mC.Parse.test.txt"
print("[PROGRESS] The output saved to %s" % outf)
df2.to_csv(outf, sep="\t", index=False)

"""
fcnt = 0
lib = {}
mxs = []
for inp in inp_ls:
  inp = inp.strip()
  head = inp.split(".5mC")[0]
  lib[head] = []

  ccnt = 0
  for i in open(inp):
    i = i.strip().split()
    ich = i[0]
    ict = int(i[1])
    if ich == ch:
      if st <= ict <= ed:
        mC = int(i[3])
        NmC = int(i[4])
        try:
          rate = round(mC/(NmC+mC),4)
        except TypeError:
          print(inp)
          print(i)
          sys.exit()
        lib[head].append(rate)
        ccnt += 1
      else: pass
    else: pass
  fcnt += 1
  print("[PROGRESS] %s/%s files processed, eluted %s C-sites" % (fcnt, len(inp_ls), ccnt))
  mxs.append(ccnt)

keys = sorted(list(lib.keys()))
mxs = max(mxs)


################ Add "NA" cells to fit to a matrix
for k in keys:
  origln = len(lib[k])
  needln = mxs - origln
  na = []
  for a in range(needln):
    na += ["NA"]
  lib[k] = lib[k] + na
  #print("[PROGRESS] Length of %s: %s" % (k, len(lib[k])))

############## Generate Output file using Pandas
out = pandas.DataFrame.from_dict(lib, orient='columns')
outf = "5mC.Parse.tmp.txt"
print("[PROGRESS] The output saved at %s" % outf)
out.to_csv(outf, sep="\t", index=False)
print()
"""
