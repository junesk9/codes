#!/usr/bin/env python3


"""
Count SNPs of bin (default 10k) from a VCF file
only "PASS" annotated SNPs will be counted.

[USAGE]
VCF.SNP-count.py [-s|] [VCF]

[OPTION]
-s : Turn-on the count-by-sample mode

[CHANGE LOG]
15.12.31 ADD count-by-sample mode, flag "-s" to turn on the mode
"""


import sys



################################ Arguments
vcf = sys.argv[-1]
abin = 1000000

if "-s" in sys.argv:
	sample_mode = True
else:
	sample_mode = False


############################### Genome distribution mode
if sample_mode == False:
	tmp = {}
	for i in open(vcf):
		i = i.strip().split()
		if i[0].startswith("#"):
			pass
		#elif i[6] == "PASS":
		else:
			ch   = i[0]
			site = int(i[1])
			try:
				tmp[ch].append(site)
			except:
				tmp[ch] = [site]
	vcf = tmp
	count = {}
	for i in vcf.keys():
		count[i] = []
		maxv = int(max(vcf[i])/abin) + 1
		for l in range(maxv):
			count[i].append(0)
	for i in count.keys():
		for site in vcf[i]:
			idx = int(site/abin)
			count[i][idx] += 1

	maxl = max([len(i) for i in count.values()])

	for l in count.keys():
		if len(count[l]) < maxl:
			for i in range(maxl-len(count[l])):
				count[l].append("NA")
		else:
			pass
	
	keys = sorted(count.keys())
	print(*keys,sep="\t")
	for l in range(maxl):
		for key in keys:
			print(count[key][l],end="\t")
		print()

######################################################### per Sample counding mode
elif sample_mode == True:

	count = {}
	sample_id = {}
	for i in open(vcf):
		i = i.strip()
		if i.startswith("##"):
			pass
		elif i.startswith("#C"):
			i = i.split()[9:]
			for en,sample in enumerate(i):
				count[en] = {}
				sample_id[en] = sample
		else:
			i = i.split()[9:]
			for en, snp in enumerate(i):
				snp = snp.split(":")[0].split("/")
				if snp[0] == ".":
					snp = "m"
				else:
					snp = int(snp[0])+int(snp[1])
				try: # Make sure to declare the count[en] as "dict" class
					foo = type(count[en])
				except NameError:
					count[en] = {}
				#end of try
				try:
					count[en][snp] += 1
				except:
					count[en][snp] = 1
				#end of try

	header = []
	for i in count.values():
		#print(i.keys())
		#sys.exit()
		i = list(i.keys())
		header = header + i
	header = list(set(header))
	header2 = ["#sample"] + header
	print(*header2, sep="\t")
	for en in sorted(count.keys()):
		print(sample_id[en], end="\t")
		for i in header:
			try:
				print(count[en][i],end="\t")
			except KeyError:
				print("NA",end="\t")
		print()
			
