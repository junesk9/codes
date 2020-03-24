#!/usr/bin/env python3

"""
SNP ratio extractor from bulked F2 gDNA seq data
based on a GATK-analyzed multi-sample VCF file.



junesk9 08th.Jun 2016
===Change log====
FEB13th 2017	Add FDR calculation; more compitability for other samples
OCT28th 2016	Add range options for comparison, with [st_peak] and [ed_peak] 
           	if range set to 25-50, the sum of two columns T25-T50 : B25-B50 shall be calculated.
"""


import sys, time
import scipy.stats as ss
import numpy as np
import pandas as pd # ver 0.18.1
import statsmodels.sandbox.stats.multicomp as mc # ver 0.8.0

ctime = time.time()

############################ Methods
def count(st):
	try:
		tmp = st.split(":")[1].split(",")
		tmp = [int(tmp[0]),int(tmp[1])]
	except:
		print(st)
		sys.exit()
	
	return tmp
def ratio(ls):
	try:
		tmp = ls[0]/sum(ls)
	except ZeroDivisionError:
		tmp = 0.5
	tmp = round(tmp,3)
	#tmp = str(tmp)

	return tmp

def Cumulative_sum(line,ls):
	ans = [0,0]
	for idx in ls:
		tmp = count(line[idx])
		ans = [ans[0]+tmp[0],ans[1]+tmp[1]]
	return ans


############## Threshold
diff_th = 0.3
dp_th = 8
pv_th = 0.1

st_peak = "50"
ed_peak = "75"

slide_mode = False
tri_mode = True ##Meaning triangular comparison?
############


############################# Arguments
x,y,z = 0,0,0
idx = {"NEUT":19,"T25":9,"T50":12,"T75":13,"T100":14,"B25":15,"B50":16,"B75":17,"B100":18,"P1":10,"P2":11}
mvcf = sys.argv[-1]
if len(sys.argv) == 1:
	try:
		mvcf = open("KMO.merged.SNP.raw.vcf")
		mvcf = "KMO.merged.SNP.raw.vcf"
	except:
		print("[ERROR] No input VCF file found! try again")
		print("__DOC__")
		sys.exit()


tp_st_idx = idx["T"+st_peak]
bt_st_idx = idx["B"+st_peak]
tp_ed_idx = idx["T"+ed_peak]
bt_ed_idx = idx["B"+ed_peak]
#print(tp_st_idx,tp_ed_idx,bt_st_idx,bt_ed_idx)

"""
for i in sorted(idx.keys()):
	cond = i
	cont_idx = idx[i]
"""

#tri_mode = True
if tri_mode == True:
	#top = "T100"
	#btm = "B" + top[1:]
	#tp_idx = idx[top]
	#bt_idx = idx[btm]
	if st_peak == ed_peak:
		top = "T"+st_peak
		btm = "B"+st_peak
	else:
		top = "T"+st_peak+"-"+ed_peak
		btm = "B"+st_peak+"-"+ed_peak
	#end of if
	n_idx = [19]

	#tp_idx = list(range(15,tp_idx+1)) ## Count All genotypes from T25; for ex, given T100 to count T25 to T100
	#bt_idx = list(range(11,bt_idx+1)) ## Count All genotypes from B25; for ex, given B100 to count B25 to B100
	tp_idx = list(range(tp_st_idx,tp_ed_idx+1))
	bt_idx = list(range(bt_st_idx,bt_ed_idx+1))
	#print(tp_idx, bt_idx)

	c = str(diff_th)[-1]
	p = str(pv_th)[-1]
	prefix = "BayRRS.SNP-FRAC.DIFF.%s.d%sc%sp%s.%s.txt" % (top, dp_th, c, p, str(time.time())[-3:])
	#out = open(prefix,"w")
	#print("SNP","CHR","POS","Bay0","RRS7",btm,"NEUT",top,"DIFF","PVAL","FDR",sep="\t",file=out)
	print("[PROGRESS] Start calculate SNP freq difference, output written to %s" % prefix)
	columns = ["SNP","CHR","POS","Bay0","RRS7",btm,"NEUT",top,"DIFF","PVAL"]
	df = pd.DataFrame(columns=columns)
	for i in open(mvcf):
		if not i.startswith("#"):
			z += 1
			i = i.strip().split()
			ch = i[0]
			pos = int(i[1])
			snp = ch+":"+str(pos)+":"+i[3]+":"+i[4]
			cvi = i[idx["P1"]].split(":")[0]
			go7 = i[idx["P2"]].split(":")[0]
			neut = i[idx["NEUT"]].split(":")[0]
			
			if cvi != go7 and "." not in cvi+go7+neut:
				cvi = count(i[idx["P1"]])
				go7 = count(i[idx["P2"]])
				neut = count(i[idx["NEUT"]])
				#top = count(i[tp_idx[-1]])
				#btm = count(i[bt_idx[-1]])

				top = Cumulative_sum(i,tp_idx)
				btm = Cumulative_sum(i,bt_idx)
				"""
				if pos == 32815:
					print(top,btm)
				"""
				diff = abs(ratio(top)-ratio(btm))
				pval = ss.fisher_exact([top,btm])[-1]
				ratio_neut = (ratio(top)-ratio(neut)) * (ratio(btm)-ratio(neut)) 
				#the ratio_neut should be negative when "neut" positioned at middle of top and bottom values.
				if diff >= diff_th and sum(neut) >= dp_th and sum(top) >= dp_th and sum(btm) >= dp_th and pval <= pv_th and ratio_neut < 0:
					if ratio(cvi) == 0 and ratio(go7) == 1:
						x += 1
						pd_line = pd.DataFrame([[snp,ch,pos,ratio(cvi),ratio(go7),ratio(btm),ratio(neut),ratio(top),diff,pval]],columns=columns)
						df = df.append(pd_line, ignore_index=True)
						#print(snp,ch,pos,ratio(cvi),ratio(go7),ratio(btm),ratio(neut),ratio(top),diff,pval,sep="\t",file=out)
					elif ratio(cvi) == 1 and ratio(go7) == 0: ## Convert all Cvi & Got7 geno.freq value to 0 & 1, respectively
						x += 1
						cvi = ratio(cvi) - 1 #convert cvi_ratio 1 to 0
						go7 = ratio(go7) + 1 # convert go7_ratio 0 to 1
						neut = round(1 - ratio(neut),3) ## convert 0.8 to 0.2
						top = round(1 - ratio(top),3)
						btm = round(1 - ratio(btm),3)
						#print(snp,ch,pos,cvi,go7,btm,neut,top,diff,pval,sep="\t",file=out)
						pd_line = pd.DataFrame([[snp,ch,pos,cvi,go7,btm,neut,top,diff,pval]], columns=columns)
						df = df.append(pd_line, ignore_index=True)
					else: pass

			if z == 100000:
				y += z
				z = 0
				ttime = time.time()-ctime
				print("[PROGRESS] %s %s lines processed, %s SNPs retained, takes %sm %ss" % (prefix, y,x,int(ttime/60),int(ttime%60)))

	print("[PROGRESS] Tottaly %s lines processed, Start FDR calculate" % y)
	pval = np.array(df.PVAL.tolist())
	fdr = mc.multipletests(pval, method="fdr_bh")
	df["FDR"] =fdr[1]
	df.to_csv(prefix, sep="\t", index=False)

	ttime = time.time()-ctime
	#out.close()
	print("[PROCESS] %s file generated" % prefix)
	print("[DONE] Totally %s lines processed, %s SNPs retained, takes %sm %ss" % ((y+z),x,int(ttime/60),int(ttime%60)))
	sys.exit()			
			
				

"""
####################################
##########

cond = "B100"
c_idx = idx[cond]
ctrl = "NEUT"
n_idx = [idx[ctrl]]

if cond.startswith("T"):
	c_idx = list(range(15,c_idx+1))
elif cond.startswith("B"):
	c_idx = list(range(11,c_idx+1))
if ctrl.startswith("T"):
	n_idx = list(range(15,n_idx[0]+1))
elif ctrl.startswith("B"):
	n_idx = list(range(11,n_idx[0]+1))


## If you want to compare particular lines
#n_idx = [14]
#c_idx = [18]

c = str(diff_th)[-1]
p = str(pv_th)[-1]
prefix = "KMO.SNP.FREQ.%s.d%sc%sp%s.txt" % (cond, dp_th, c,p)
out = open(prefix,"w")
#out = open("t100-b100.diff.txt","w")
print("CHR","POS","Cvi0","Go7","NEUT",cond,"DIFF","PVAL",sep="\t",file=out)
print("[PROGRESS] Start calculate SNP frequency, output written to %s" % prefix)


chrs, ch_ln = [],[]

for i in open(mvcf):
	if not i.startswith("#"):
		z += 1
		i = i.strip().split()
		ch = i[0]
		pos = int(i[1])
		cvi = i[9].split(":")[0] 
		go7 = i[10].split(":")[0]
		
		##assess chromosomes number and length
		if ch_ln == chrs == []:
			chrs.append(ch)
			ch_ln.append(pos)
		elif ch == chrs[-1]:
			ch_ln[-1] = pos if ch_ln[-1]<pos else ch_ln[-1]
		elif ch != chrs[-1]:
			chrs.append(ch)
			ch_ln.append(pos)


		if cvi != go7 and "." not in cvi+go7: ## Both are not missing & different to each other; "." from "./."
			cvi = count(i[9])
			go7 = count(i[10])
			neut = Cumulative_sum(i,n_idx)
			t25 = Cumulative_sum(i,c_idx)
			
			#neut = [0,0]
			#for idx in n_idx: # Cumulative ratio calculation
				#tmp = count(i[idx])
				#neut = [neut[0]+tmp[0],neut[1]+neut[1]]
			#end of for
			#t25 = [0,0]
			#for idx in c_idx: ## Cumulative ratio calculation
				#tmp = count(i[idx])
				#t25 = [t25[0]+tmp[0],t25[1]+tmp[1]]
			#end of for
			diff = abs(ratio(neut)-ratio(t25))
			pval = ss.fisher_exact([neut,t25])[-1]
			if diff >= diff_th and sum(neut) >= dp_th and sum(t25) >= dp_th and pval <= pv_th: #difference cut-off
				x += 1
				#print(diff)
				if ratio(cvi) == 0 and ratio(go7) == 1:
					print(ch,pos,ratio(cvi),ratio(go7),ratio(neut),ratio(t25),diff,pval,sep="\t",file=out)
				elif ratio(go7) == 0 and ratio(cvi) == 1:
					cvi = ratio(cvi) - 1 #convert cvi_ratio 1 to 0
					go7 = ratio(go7) + 1 # convert go7_ratio 0 to 1
					neut = round(1 - ratio(neut),3) ## 0.8 to 0.2
					t25 = round(1 - ratio(t25),3)
					diff = 1-diff
					print(ch,pos,cvi,go7,neut,t25,diff,pval,sep="\t",file=out)
		else: pass
		#end of if
		if z == 100000:
			y += z
			z = 0
			ttime = time.time()-ctime
			print("[PROGRESS] %s %s lines processed, %s SNPs retained, takes %sm %ss" % (prefix, y,x,int(ttime/60),int(ttime%60)))
		else: pass
		#end of if
ttime = time.time()-ctime
out.close()
print("[PROCESS] %s file generated" % prefix)
print("[DONE] Totally %s lines processed, %s SNPs retained, takes %sm %ss" % ((y+z),x,int(ttime/60),int(ttime%60)))

"""


######Further SNP-counting by slides steps

if slide_mode == True:
	window = 1 ##megabasepairs
	slide = 10 ##kilobasepairs


	print()
	print("[PROGRESS] Start window sliding analysis of %s, window size: %sM, slide gap: %sK" % (prefix,window,slide))
	ch_ln = [int(i/(slide*1000))+1 for i in ch_ln]
	print("[PROGRESS] Reference has %s chr, %s windows to be analyzed" % (len(chrs), sum(ch_ln)))

	prefix2 = ".".join(prefix.split("."))[:-2]+".slide.txt"
	fi = [i.strip().split() for i in open(prefix)]
	header = fi[0]
	fi = fi[1:] # w/o header
	fi = pd.DataFrame(fi, columns=header)
	fi.POS = fi.POS.astype(float).fillna(0.0) ## Covert Pandas data to numeric
	fi.PVAL = fi.PVAL.astype(float).fillna(1.0)


	pval =["e2","e3","e4","e5","e6"]
	pval_dic = {"e2":0.01,"e3":0.001,"e4":0.0001,"e5":0.00001,"e6":0.000001}

	out2 = open(prefix2,"w")
	print("SNP","CHR","POS","e2","e3","e4","e5","e6",sep="\t",file=out2)
	x,y,z = 0,0,0
	for en,ch in enumerate(chrs):
		for win in range(ch_ln[en]):
			x += 1
			mv = win*slide*1000
			r_min, r_max = 0+mv, (window*1000000)+mv
			print(snp,ch,r_min,sep="\t",end="\t",file=out2)
			for i in pval:
				ev = pval_dic[i]
				tmp = fi.query("CHR == @ch & PVAL < @ev & POS < @r_max & POS >= @r_min")
				#print(tmp)
				tmp = tmp.iloc[:,5].tolist() ## extact 6th column (t25) to list
				tmp = [float(i) for i in tmp] ## convert the data to float
				#print(tmp)
				tmp = np.mean(tmp) if len(tmp)>9 else "nan"
				#print(tmp)
				#pval_dic[i].append(tmp)
				tmp = "0" if tmp == "nan" else tmp
				print(tmp,end="\t",file=out2)
			print(file=out2)
		print("[PROGRESS] Sliding %s chr done, %s windows evaluated" % (ch, ch_ln[en]))
			
	out2.close()

