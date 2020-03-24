#!/usr/bin/env python3

"""
QTL-seq like analyzing with a globla VCF file

Fisrt filtered the input VCF file, just check the parental SNP rates in two parents are distinctive (0,1)

and in the second step calculate SNP ratio in a sliding window (default: 1M by 10k)
the mean of SNP-difference between top-group and btm-group assessed,
with multiple significancy values (1% 5% 95% 99%) are also given. (by sum(top+btm) - btm) 


now only singl input (a gVCF file from GATK)

Apr 4th 2017 Junesk9

"""

import sys, time
import scipy.stats as ss
import numpy as np
import pandas as pd # ver 0.18.1
import statsmodels.sandbox.stats.multicomp as mc # ver 0.8.0

ctime = time.time()

################### Arguments
x,y,z = 0,0,0
#idx = {"NEUT":19,"T25":9,"T50":12,"T75":13,"T100":14,"B25":15,"B50":16,"B75":17,"B100":18,"P1":10,"P2":11}
idx = {"NEUT":19,"T25":15,"T50":16,"T75":17,"T100":18,"B25":11,"B50":12,"B75":13,"B100":14,"P1":9,"P2":10}
args = sys.argv
mvcf = args[-1]
#mvcf = "BR-F2.2nd.FILT.NAMED2.recode.vcf" ## temp
mvcf = "KMO.merged.filt2.vcf.recode.vcf"

diff_th = 0.3
dp_th = 8
pv_th = 0.1

st_peak = "50"
ed_peak = "50"

window = 1 ## megabasepairs
slide = 10 ## kilobasepairs

################################ Run


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

def Cumul_sum(doubled_list):
	ans = [0,0]
	for m in doubled_list:
		m = [int(m[0]),int(m[1])]
		ans = [ans[0]+m[0],ans[1]+m[1]]
	return ans

def ListedDelta(top, btm):
	t, b = [], []
	
	for l in top:
		tmp = ratio(l)
		t.append(tmp)
	for l in btm:
		tmp = ratio(l)
		b.append(tmp)
	ans = []
	for en,i in enumerate(t):
		tmp = float(t[en]-b[en])
		ans.append(tmp)
	
	return sorted(ans)

def makeHist():
        arg = {"g":""}
        if "-g" in args:
                gidx = args.index["-g"]+1
                arg["g"] = args[gidx]

                gdic = {}
                gn = arg["g"]
                for i in open(gn):
                        i = i.strip()
                        if i.startswith(">"):
                                ch = i.split()[0][1:]
                                gdic[ch] = []
                        else:
                                gdic[ch] + [i]
                for i in gdic.keys():
                        gdic[i] = "".join(gdic[i])
                        gdic[i] = len(gdic[i])
        else:
                gdic = {"2":19698289,
                        "1": 30427671,
                        "4": 18585056,
                        "5": 26975502,
                        "3": 23459830} #TAIR10 genome 

        hist = {}
        for k in gdic.keys():
                win_size = 1+int(int(gdic[k])/(slide*1000))
                hist[k] = win_size		
        return hist

#print(makeHist())
#sys.exit()


def Pval(ls):
	cnt = len(ls)
	
	if not cnt == 0:
		ls = sorted(ls)
		perc1 = ls[int(cnt*0.01)]
		perc5 = ls[int(cnt*0.05)]
		perc95 = ls[int(cnt*0.95)]
		perc99 = ls[int(cnt*0.99)]
	else:
		perc1,perc5,perc95,perc99 = "NA","NA","NA","NA"

	return [perc1,perc5,perc95,perc99]



#############################################################


################################ Run


tp_st_idx = idx["T"+st_peak]
bt_st_idx = idx["B"+st_peak]
tp_ed_idx = idx["T"+ed_peak]
bt_ed_idx = idx["B"+ed_peak]
n_idx = idx["NEUT"]

if st_peak == ed_peak:
	top = "T"+st_peak
	btm = "B"+st_peak
else:
	top = "T"+st_peak+"-"+ed_peak
	btm = "B"+st_peak+"-"+ed_peak


tp_idx = list(range(tp_st_idx,tp_ed_idx+1))
bt_idx = list(range(bt_st_idx,bt_ed_idx+1))
#print(tp_idx, bt_idx)

c = str(diff_th)[-1]
p = str(pv_th)[-1]
tm = str(time.time())[-3:]
#prefix = "BayRRS.QTLseq.res1.%s.d%sc%sp%s.%s.txt" % (top, dp_th, c, p, tm)
prefix2 = "BayRRS.QTLseq.%sd%s.win%sM.sld%sK.txt" % (top, dp_th, window, slide)
print("[PROGRESS] Start QTLseq analysis of %s" % mvcf)
#columns = ["CHR","Win_st","Win_ed","SNP_count","Mean_P1","Mean_P2","Mean_delta","1perc","5perc","95perc","99perc"]
columns = ["CHR","POS","P1","P2","NEUT",top,btm]
#df = pd.DataFrame(columns=columns)
listed = []
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

		if "." not in cvi+go7+neut and go7 != cvi:
			cvi = count(i[idx["P1"]])
			go7 = count(i[idx["P2"]])
			neut = count(i[idx["NEUT"]])
			top = Cumulative_sum(i,tp_idx)
			btm = Cumulative_sum(i,bt_idx)
			#pv_top = ss.fisher_exact([top,neut])[-1]
			#pv_btm = ss.fisher_exact([btm,neut])[-1]			
			
			if sum(top) >= dp_th and sum(btm) >= dp_th:
				if ratio(cvi) == 0 and ratio(go7) == 1:
					lst = [ch,pos,cvi,go7,neut,top,btm]
					listed += [lst]
					x += 1
				elif ratio(cvi) == 1 and ratio(go7) == 0: ## Set value of cvi(weak accs) to 0, value of go7(strong accs) to 1 
					cvi = [cvi[1],cvi[0]]
					go7 = [go7[1],go7[0]]
					neut = [neut[1],neut[0]]
					top = [top[1],top[0]]
					btm = [btm[1],btm[0]]
					lst = [ch,pos,cvi,go7,neut,top,btm]
					listed += [lst]
					x += 1
				else: pass
			"""
			if ratio(cvi) in [0,1] and ratio(go7) in [0,1] and sum(top) >= dp_th and sum(btm) >= dp_th:  #pv_top < pv_th and pv_btm < pv_th:
				lst = [ch,pos,cvi,go7,neut,top,btm]
				listed += [lst]
				#pd_line = pd.DataFrame([[ch,pos,cvi,go7,neut,top,btm]], columns=columns)
				#df = df.append(pd_line, ignore_index=True)
				x += 1
			"""
	if z == 500000:
		y += z
		z = 0
		ttime = time.time()-ctime
		print("[PROGRESS] %s lines processed, %s retained, takes %sm %ss" % (y, x,int(ttime/60),int(ttime%60)))
print("[PROGRESS] Tottaly %s lines processed, converting to pandas.DataFrame" % (y+z))
df = pd.DataFrame(listed, columns=columns)
#print("[PROGRESS] SNP collection done, the intermediates written to %s" % prefix) 
#df.to_csv(prefix, sep="\t", index=False)

#print(df)
#sys.exit()

#################
print("[PROGRESS] Start %s-window %s-sliding" % (str(window)+"M", str(slide)+"K"))
#prefix2 = "BayRRS.QTLseq.%sd%s.win%sM.sld%sK.txt" % (top, dp_th, window, slide)
columns = ["CHR","Win_st","Win_ed","SNP_count","Mean_P1","Mean_P2","Mean_delta","1perc","5perc","95perc","99perc"]
#columns = ["CHR","Win_st","Win_ed","SNP_count","Mean_P1","Mean_P2","Mean_delta","pval_neut","pval_sum"]
#df2 = pd.DataFrame(columns=columns)

gdic = makeHist()

z = 0
listed = []
chrs = list(sorted(gdic.keys()))
for k in chrs:
	print("[PROGRESS] Start CHR %s by %s slidings" % (k, gdic[k]))
	wnd = window*1000000
	w,x,y = 0,0,0
	for r in range(gdic[k]):
		r = r * slide * 1000
		rng = [r, r+wnd]
		#print(rng)
		st = rng[0]
		ed = rng[1]

		subset = df.query("CHR == @k & POS >= @st & POS < @ed")
		#print(subset)
		#sys.exit()
		p1 = subset.iloc[:,2].tolist()
		p2 = subset.iloc[:,3].tolist()
		neut = subset.iloc[:,4].tolist()
		top = subset.iloc[:,5].tolist()
		btm = subset.iloc[:,6].tolist()
		tot = []
		for en,i in enumerate(btm):
			t = top[en]
			tmp = [i[0]+t[0], i[1]+t[1]]
			tot.append(tmp)


		#rat_p1 = [ratio(i) for i in p1]
		#rat_p2 = [ratio(i) for i in p2]
		#rat_top = [ratio(i) for i in top]
		#rat_btm = [ratio(i) for i in btm]
		#rat_neut = [ratio(i) for i in neut]
		
		ld = ListedDelta(top, btm)
		ld_neut = ListedDelta(neut, btm)
		ld_sum = ListedDelta(tot, btm)

		mean_p1 = ratio(Cumul_sum(p1))
		mean_p2 = ratio(Cumul_sum(p2))
		delta = np.mean(ld) if len(ld) != 0 else "NA"
		pval_neut = Pval(ld_neut)
		pval_sum = Pval(ld_sum)

		#print(mean_ld, pval_neut, pval_sum)	
		#sys.exit()	

		#print("ld: ", ld[:5])
		#print("ld_neut: ", ld_neut[:5])
		#print("ld_sum: ", ld_sum[:5])
		#sys.exit()
		
		"""
		rat_neut = sorted(rat_neut)
		perc1 = rat_neut[int(len(rat_neut)*0.01)]
		perc5 = rat_neut[int(len(rat_neut)*0.05)]
		perc95 = rat_neut[int(len(rat_neut)*0.95)]
		perc99 = rat_neut[int(len(rat_neut)*0.99)]
		print(perc1, perc5, perc95, perc99)
		#sys.exit()		

		mean_p1 = ratio(Cumul_sum(p1))
		mean_p2 = ratio(Cumul_sum(p2))
		mean_top = ratio(Cumul_sum(top))
		mean_btm = ratio(Cumul_sum(btm))
		mean_neut = ratio(Cumul_sum(neut))
		delta = float(mean_top - mean_btm)

		print(mean_p1)
		print(mean_p2)
		print(mean_top)
		print(mean_btm)
		print(mean_neut)
		

		
		neut_ary = np.array([ratio(i) for i in neut])
		top_ary = np.array([ratio(i) for i in top])
		btm_ary = np.array([ratio(i) for i in btm])
		#dlt_ary = np.array([abs(i) for i in ListedDelta(top, btm)])
		pval_top = ss.wilcoxon(neut_ary, top_ary)
		pval_btm = ss.wilcoxon(neut_ary, btm_ary)
		print("neut_array: ",len(neut_ary),neut_ary[:5])
		#print("dlt_array: ",len(dlt_ary), dlt_ary[:5])
		print("pval_neut: ",pval_top, pval_btm)
		
		
		#pval_neut = ss.fisher_exact([     ])[-1] 
		#pval_sum = 

		lb = ListedDelta(top, btm)
		print(np.median(lb))
		print(np.mean(lb))
		print(delta)
		sys.exit()
		if lb != []:
			p_ls = []
			for en,n in enumerate(neut):
				if delta > 0:
					pval = ss.fisher_exact([top[en],neut[en]])[-1]
				elif delta < 0:
					pval = ss.fisher_exact([btm[en],neut[en]])[-1]
				else:
					pval = 1
				p_ls.append(pval)
			p_ls = sorted(p_ls)
			perc1 = p_ls[int(len(p_ls)*0.01)]
			perc5 = p_ls[int(len(p_ls)*0.05)]
			perc95 = p_ls[int(len(p_ls)*0.95)]
			perc99 = p_ls[int(len(p_ls)*0.99)]

		else:
			w += 1
			perc1, perc5, perc95, perc99 = "NA","NA","NA","NA"
		"""
		perc1, perc5, perc95, perc99 = pval_sum[0], pval_sum[1], pval_sum[2], pval_sum[3]
		pd_line = [k,st,ed,len(p1),mean_p1,mean_p2,delta,perc1,perc5,perc95,perc99]
		listed.append(pd_line)	
		#pd_line = pd.DataFrame([[ch,st,ed,len(p1),mean_p1,mean_p2,delta,perc1,perc5,perc95,perc99]], columns=columns)
		#df2 = df.append(pd_line, ignore_index=True)
		x += 1
		if x == 1000:
			y += x
			z += x
			x = 0
			ttime = time.time()-ctime
			print("[PROGRESS] CHR %s, %s lines processed, takes %sm %ss" % (k, y,int(ttime/60),int(ttime%60)))
	print("[PROGRESS] CHR %s done, %s line(s) without value" % (k,w))
print("[PROGRESS] All CHRs processed, tottally %s lines " % (z+x))
df2 = pd.DataFrame(listed, columns=columns)
print("[PROGRESS] The result saved into %s file" % prefix2)
print(df2)
df2.to_csv(prefix2, sep="\t", index=False)

				


