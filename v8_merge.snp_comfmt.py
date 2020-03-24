#!/usr/bin/env python3

"""
As the second step from the resultant "*.txt.2" of a.var2comfmt.py;
Merging common-formatted snp data to single file "snp.merged.txt"
Equip the ref sequence validation steps, matching nt at the snp data vs. nt indicated by indexing on the tair10.fa


[USAGE]
-at the folder containing snp data files & tair10.genome.fa
>python3 b.merge.comfmt.py -ref [reference genome fasta]

--The common snp format file should have name as ACCS.comFmt.txt.2 & contain 4 columns with tab-deliminitor
[chr][site][ref][alt]

--The output file contains 5+[#accs] columns
[chr][site][snp_id][ref][alt][accs1][accs2][accs3]....

[alt] indicates all unique alternative (different to the reference) nt from all merging snp data
      "*" refers single base deletion
      "-" refers no data; whether missing or no-snp is unknown

junesk9

====LOG====
Nov12th2015	costum reference genome available
Oct7th2015 	First build
"""
import sys, glob

############################### Arguments

txt = glob.glob("*.txt.2")
txt = sorted(list(txt))
vrs = [i.split(".")[0] for i in txt]

args = sys.argv[1:]
try:
	ref_idx = args.index("-ref")+1
	ref = args[ref_idx]
except:
	print(__doc__)
	sys.exit(1)

print("""
	Merging multi-SNP files 

	# of SNP data: %s
	reference    : %s

	""" % (len(txt),ref))


################################## Prepare reference TAIR10 genome for following SNP validation steps.
ref = open(ref)
fa = {}
for i in ref:
	i = i.strip().split()[0]
	if i.startswith(">"):
		name = i[1:]
		fa[name] = []
	else:
		fa[name].append(i)
for i in fa.keys():
	fa[i] = "".join(fa[i])
	


########################### Merging SNP sites
sites = []

for fi in txt:
	for i in open(fi):
		i = i.strip().split("\t")
		site = i[0]+":"+i[1]
		sites.append(site)
sites = sorted(list(set(sites)))
print("Total collected SNPs are %s lines" % len(sites))

site = []
for i in sites:
	site.append([i])
sites = site


############################## Make list dic [var1{ch:site : ref:alt, ...}, var2]
ls = {}
for fi in txt:
	name = fi.split(".")[0]
	ls[name] = {}
	for i in open(fi):
		i = i.strip().split("\t")
		site = i[0]+":"+i[1]
		if i[3] == "-":
			i[3] = "*"
		seq = i[2]+":"+i[3]
		if i[3] == "N" or i[2] == "N":
			seq = "-"
		ls[name][site] = seq

x,y = 0,1 ## Counting parts
for st in sites:
	for vr in vrs:
		if st[0] in ls[vr].keys():
			st.append(ls[vr][st[0]])
		else:
			st.append("-")
		x += 1
		if x == 10000000: # Counting parts
			print("%s lines processed" % (y*10000000))
			x = 0
			y += 1

print("%s lines processed" % (x+(y-1)*10000000) )
#print(*sites[:10],sep="\n")
### Contents of "sites"
#['Chr1:10000121', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'T:C', '-', '-']
#['Chr1:10000153', '-', '-', '-', '-', '-', '-', 'C:T', 'C:T', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']
#['Chr1:10000173', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'A:C', '-', '-', 'A:C', '-', 'A:C', '-', '-', '-', '-', '-']
tmp_out = open("merged.NotFiltered.txt","w")
for i in sites:
	print(*i,sep="\t",file=tmp_out)
tmp_out.close()

print("Generate merged.NotFiltered.txt Done")

############################# Parse the list to more convinient view & validate the SNP concordance

sites = [i.strip().split() for i in open("merged.NotFiltered.txt")]
print("Load Not filtered data done")

snp = []

for j in sites:
	tmp = j[0].split(":")
	seqs = j[1:]
	ref,alt = [],[]
	for k in seqs:
		if ":" in k:
			seq = k.split(":")
			ref.append(seq[0])
			alt.append(seq[-1])
	ref = ",".join(sorted(set(ref)))
	alt = ",".join(sorted(set(alt)))
	tmp = tmp + [ref,alt] + seqs
	snp.append(tmp)
			
#print(*snp[:10],sep="\n")

er = 0
passed,notPassed = [],[]
for i in snp:
	ch = i[0][-1]
	idx = int(i[1])-1
	ref = i[2]
	alt = i[3]
	try:
		if len(ref) == 1 and ref == fa[ch][idx]:
			passed.append(i)
		else: 
			notPassed.append(i)
	except IndexError:
		er += 1
		#print(i)
		pass
	except KeyError:
		er += 1
		pass
print(len(passed),len(notPassed),er)

snp_count = []
for i in vrs:
	snp_count.append(0)
#snp_count = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

out = open("snp.merged.filtered.txt","w")
header = ["#CHR","BP","SNP_ID","REF","ALT"] + vrs
print(*header,sep="\t",file=out)
for i in passed:
	info = i[:4]
	snp_id = ":".join(info)
	print(info[0],info[1],snp_id,info[2],info[3],sep="\t",end="\t",file=out)
	i = i[4:]
	for en,j in enumerate(i):
		if j != "-":
			snp_count[en] += 1
		print(j[-1],end="\t",file=out)
	print(file=out)
print(snp_count)
out.close()

