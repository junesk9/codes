#!/usr/bin/env python3
"""
SNP data merger

From multiple snp.txt files (modified vcf file; the output of vcf.parser.py, looks like)
make single snp file contains the every allele info for each snp site
equips SNP_quality score filteration, adopt from .VCF file column 6

## Input file format
[gene] [site] [chr] [gene_st] [gene_ed] [ref] [alt] [alt_count] [ref_count] [freq] [snp_qual]

## Ouput format
[gene] [site] [chr] [gene_st] [gene_ed] [ref] [alt(sample1)] [alt(smaple2)] ...

If the snp_qual was relatively low or the snp info was absent in the particular sample file,
the [alt] appears as [ref] with lower case.

If the SNP frequency was lower than 0.95, recognized the SNP as heterozygous,
the [alt] appears as [alt/ref] form.


ToshiMori9 May. 31st, 2015
"""
import sys, glob, pickle, time


cutoff = 60
prefix = "out.qual_"+str(cutoff)+".snp"
ls = sorted(glob.glob("[6-9]*txt"))

site_ls = []
for i in ls:
	fi = open(i).readlines()[1:]
	for j in fi:
		j = j.strip().split("\t")
		if float(j[-1]) >= cutoff: ## SNP_qual filtering
			site = ":".join(j[0:6])
			site_ls.append(site)

site = sorted(list(set(site_ls)))
print("SNP sites collected, totally %s sites" % len(site))
print("%s VCF files will be parsed, %s" % (len(ls), time.ctime()))


## Make dictionary
dic = {}
for i in site:
	key = ":".join(i.split(":")[0:5])
	dic[key] = [i.split(":")[-1]]


header = ["gene","site","chr","gene_st","gene_ed","ref"]
for en,fi in enumerate(ls):
	fis = open(fi).readlines()[1:]

	## Output hedaer formation
	head = fi.split(".")[0]
	header.append(head)

	x, y = 0, 1 ## Counting body
	for i in fis:
		i = i.strip().split("\t")
		loci = ":".join(i[0:5])
		ref = i[5].lower()
		alt = i[6]
		freq = float(i[9])
		qual = float(i[-1])
		if qual >= cutoff: #if acceptable SNP
			if freq < 0.95: # determine homo/heterozygous SNP
				alt = alt+"/"+ref
			else: pass
			try:
				dic[loci].append(alt)
			except KeyError: pass
		else:
			try:
				dic[loci].append(ref)
			except KeyError: pass
		"""
		# Counting body
		x += 1
		if x == 100000:
			print("%s, %s lines done!" % (fi, x*y))
			x = 0
			y += 1
		"""
	for k, v in dic.items():
		if not len(v) == en+2:
			if len(v) == en+1:
				dic[k].append(v[0].lower())
			else:
				print("ERROR")
				sys.exit()
	print("%s, %s lines done!, file end %s" % (fi, len(fis), time.ctime()))

"""
##Pickling
ctime = time.time()
pick = "vcf.parse"+str(ctime)+".pickle"
with open(pick, 'wb') as f:
	pickle.dump(dic, f)
"""

## Make ouput file
out = open(prefix, 'w')
for head in header:
	print(head, end="\t", file=out)
print(file=out)
for x in site:
	z = dic[":".join(x.split(":")[0:5])]
	x = x.split(":")[0:5]
	for y in x:
		print(y,end="\t", file=out)
	for a in z:
		print(a,end="\t", file=out)
	print(file=out)
print("Done! %s, %s" % (prefix, time.ctime()))
dic, site_ls = [], []
out.close()
