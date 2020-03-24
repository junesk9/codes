#!/usr/bin/env python3

#!/usr/bin/env python3

"""
2018.11.21 Junesk9

"""

import pandas as pd
import sys
import time

############  Option arguments

opts = {"-vcf":"", "-rmch":[], "-qualth":0, "-out":"", "-biallele":False, "-covth":60}
args = sys.argv[1:]
remove_ch = False
qual_th = False
outfile = False

### temp. opt
special_opt = True

if "-vcf" not in args:
	print("[ERROR] Requires at least one VCF files, try again with -vcf ")
	sys.exit()
elif "-vcf" in args:
	vidx = args.index("-vcf") + 1
	opts["-vcf"] = args[vidx]
if "-biallele" in args:
	opts["-biallele"] = True
if "-out" in args:
	outfile = True
	oidx = args.index("-out")+1
	opts["-out"] = args[oidx]
if "-rmch" in args:
	remove_ch = True
	chidx = args.index("-rmch")+1
	opts["-rmch"] = args[chidx].split(",")
if "-qualth" in args:
	qual_th = True
	qidx = args.index("-qualth")+1
	try:
		opts["-qualth"] = int(args[qidx])
	except:
		print("[ERROR] QUALITY THRESHOLD (-qualth) only accept an integer")
		sys.exit()

print("""
        ############## VCF stat analyzer ###############
        
        Input VCF     : %s
        remove chrs   : %s (%s)
        quality cut   : %s (%s)
	bialleles only: %s
	write output : %s (cut SNP qual: %s)

        ###############################################
      """ % (opts["-vcf"],remove_ch,opts["-rmch"],qual_th,opts["-qualth"], opts["-biallele"], opts["-out"], opts["-covth"]))

vcf, rmch, qth, cth = opts["-vcf"],opts["-rmch"],opts["-qualth"], opts["-covth"]
biallele, outf = opts["-biallele"], opts["-out"]
if qual_th == False and remove_ch == False:
	biallele = False
	print("[CAUTION] Since NO filter assigned, the outfile option is descreted")
else: pass

################## Load the VCF
print("[PROGRESS] Start open file %s" % vcf)
stime = int(time.time())
pdheader = ["CHR","BP","QUAL","COV"]
df, df_filt = [], []
#chs, quals, covs = [],[],[]
for v in open(vcf).readlines():
	if v.startswith("##"):
		pass
	elif v.startswith("#CHROM"):
		header = v.split()[1:]
		num_h = len(header)
	else:
		l = v.split()
		ch, bp, qual = l[0],l[1],l[5]
		genotypes = l[9:]
		genotypes = [g.split(":")[0] for g in genotypes]
		cov = sum(1 for i in genotypes if i != "./.")
		cov = round(float(cov*100/num_h),2)

		if biallele == True and "," in l[4]:
			pass
		elif special_opt == True and genotypes[0] == genotypes[1]:
			pass ## only Got != Cvi genotypes retained
		else:
			dfline = [ch,bp, qual,cov]
			df.append(dfline)
			if remove_ch == True:
				if dfline[0] in rmch:
					pass
				elif qual_th == True:
					if float(dfline[1]) < qth:
						pass
					else:
						df_filt.append(dfline)
				else:
					df_filt.append(dfline)

                        
                #df = df.append(dfline, ignore_index=True)
                #chs.append(ch)
                #quals.append(qual)
                #covs.append(cov)
                #chs = sorted(list(set(chs)))
                #quals = sorted(quals)

print(len(df), len(df_filt))
#sys.exit()

df = pd.DataFrame(df, columns=pdheader)
df_filt = pd.DataFrame(df_filt, columns=pdheader)
#print(df)
etime = int(time.time())
etime = etime - stime
print("[PROGRESS] VCF file loaded, takes %s secs" % etime)

####################### Metrics
def CountChrs(dataframe):
	chrs = dataframe['CHR'].tolist()
	uniq_ch = sorted(list(set(chrs)))

	print("Tot: %s" % len(chrs), end="\t")
	for en,c in enumerate(uniq_ch):
		ccount = sum(1 for i in chrs if i == c)
		if en+1 != len(uniq_ch):
			print("%s: %s" % (c, ccount), end="\t")
		else:
			print("%s: %s" % (c, ccount))

def CountCovs(dataframe):
	covs = dataframe['COV'].tolist()
	print("Tot: %s" % len(covs), end="\t")

	temp = covs
	for p in [10,20,30,40,50,60,70,80,90]:
		temp = [i for i in temp if i >= p]
		if p != 90:
			print(">%s: %s" % (p, len(temp)), end="\t")
		else:
			print(">%s: %s" % (p, len(temp)))
	"""
	perc90 = [i for i in covs if i >= 90]
	perc80 = [i for i in perc90 if i >= 80]
	perc70 = [i for i in perc80 if i >= 70]
	perc60 = [i for i in perc60 if i >= 60]
	perc50 = [i for i in perc50 if i >= 50]
	"""

########################## Print result

print("[RESULT] Raw stats")
print()
CountChrs(df)   
CountCovs(df)
print()
if remove_ch or qual_th == True:
	print("[RESULT] Filtered stats")
	print()
	CountChrs(df_filt)
	CountCovs(df_filt)


######################### Generate a filetered vcf file

if outfile == True and outf != "":
	df_filt = df_filt[df_filt['COV'] >= cth] ### 70 as current coverage threshold 2018NOV
	chs = df_filt['CHR'].tolist()
	bps = df_filt['BP'].tolist()
	idx = [i+":"+bps[en] for en,i in enumerate(chs)]

	print("[PROGRESS] %s lines output to %s" % (len(idx), outf))
	outf = open(outf,"w")
	xx, y = 0,0
	for i in open(vcf):
		i = i.strip()
		if i.startswith("#"):
			print(i, file=outf)
		else:
			v = i.split()
			chbp = v[0]+":"+v[1]
			if chbp in idx:
				idx.remove(chbp)
				print(i, file=outf)
				xx += 1
			else: pass
		if xx == 10000:
			xx = 0
			y += 1
			print("[PROGRESS] %s lines written" % int(y*10000))
		
	print("[DONE] Tottally %s lines written" % len(idx))
	outf.close()


"""
This code count SNP categories, by sample (in case of GATK-VCF)
SNP SNP(C-T) INDEL NON-SNP

Total count also will made.

DEC.05th.2016



import sys

vcf = sys.argv[-1]
if not vcf.endswith("vcf"):
	print("[ERROR] The input file should be single VCF file produced by GATK")
	sys.exit()
if len(sys.argv) > 2:
	print("[ERROR] Only single argument can be accepted to this code")
	sys.exit()


#x,y,z = 0,0,0
dic = {"tot":[0,0,0,0]}
for i in open(vcf):
	i = i.strip()
	if i.startswith("##"):
		pass
	elif i.startswith("#CH"):
		i = i.split()
		head = i[9:]
		for h in head:
			dic[h] = [0,0,0,0] #SNP,INDEL,C-T,Non
	else:
		i = i.split()
		#print(i[2],i[3],i[4],i[5])
		if len(i[3]) == len(i[4]) == 1:
			if i[3] == "C" and i[4] == "T":
				idx = 2
				dic["tot"][2] += 1
			elif i[3] == "G" and i[4] == "A":
				idx = 2
				dic["tot"][2] += 1
			else:
				idx = 0
				dic["tot"][0] += 1
		else:
			idx = 1
			dic["tot"][1] += 1
		i = i[9:]
		for en,s in enumerate(i):
			s = s.split(":")[0]
			if s in ["1/1","0/1","1/0"]:
				dic[head[en]][idx] += 1
			else:
				dic[head[en]][3] += 1

keys = sorted(list(dic.keys()))

for k in keys:
	print(k, end="\t")
	print(*dic[k],sep="\t")
"""
