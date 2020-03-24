#!/usr/bin/env python3

"""
Usage

python3 vC_VCF_genoExtractor [-l:10] vcf snp

[OPTIONS]
-l : range snps to collect genotypes default: 10

[LOG]
18.may17	Add whole VCF mode; more descriptive code
"""


import sys

if len(sys.argv) == 1:
	print(__doc__)
	sys.exit()

arg = {"vcf":"", "snp":"", "rng":10}
args = sys.argv[1:]
if "-vcf" in args:
	v_idx = args.index("-vcf") + 1
	arg["vcf"] = args[v_idx]
else:
	print(__doc__)
	print("[ERROR] One VCF file input is mendotary")
	sys.exit()

if "-snp" in args:
	s_idx = args.index("-snp") + 1
	arg["snp"] = args[s_idx]
	mode_slice = True
else:
	mode_slice = False

if "-rng" in args:
	r_idx = args.index("-rng") + 1
	try:
		arg["rng"] = int(args[r_idx]) if mode_slice == True else "NA"
	except:
		print("[ERROR] RANGE argument must be a integer; your input was %s" % args[r_idx])
		sys.exit()	
if "-out" in args:
	o_idx = args.index("-out")
	arg["out"] = args[o_idx]
else:
	if not arg["snp"] == "":
		arg["out"] = "-".join(snp.split(":"))+".geno.txt"
	else:
		arg["out"] = ".".join(arg["vcf"].split(".")[:-1]) + ".geno.txt"

#vcf = sys.argv[-2]
#snp = sys.argv[-1]
#prefix = "-".join(snp.split(":"))+".geno.txt"

print("""

	VCF Genotype extractor

	vcf  : %s
	Slice: %s

	snp  : %s
	range: %s (default 10)
	out  : %s

""" % (arg["vcf"], mode_slice, arg["snp"], arg["rng"], arg["out"]))

vcf, snp, rng, prefix = arg["vcf"], arg["snp"], arg["rng"], arg["out"]

#######################################################


geno = {}
header = ["#ID",]
vcf = [i.strip() for i in open(vcf)]

if mode_slice == True:
	for en,i in enumerate(vcf):
		if i.startswith("##"):
			pass
		elif i.startswith("#C"):
			pass_n = en
			i = i.split()[9:]
			for key in i:
				geno[key] = []
		else:
			i = i.split()[2]
			if snp == i:
				snp_id = en
			else: pass
else: #mode_slice == False
	for en, i in enumerate(vcf):
		if i.startswith("##"):
			pass
		elif i.startswith("#C"):
			i = i.split()[9:]
			header += i
		else:
			snp_id = i.split()[2]
			flags = i.split()[8].split(":")
			ad_idx = flags.index("AD")
			gt_idx = flags.index("GT")
			#print(flags, ad_idx, gt_idx)
			gt_perc = []
			for g in i.split()[9:]:
				g = g.split(":")
				gt = g[gt_idx]
				if gt == "./.":
					ad = "NA"
				else:	
					ad = g[ad_idx].split(",")
					#print(g, ad)
					ref = int(ad[0])
					alt = int(ad[1])
					if ref+alt == 0:
						if gt == "1/1":
							ad = "100(n)"
						elif gt == "0/0":
							ad = "0(n)"
						elif gt in ["1/0", "0/1"]:
							ad = "50(n)"
					else:
						ad = float(alt/(ref+alt)) * 100
				gt_perc.append(ad)
			geno[snp_id] = gt_perc
			
			
if mode_slice == True:
	keys = sorted(geno.keys())
	print("[PROCESS] Totally %s variants of %s individuals on calculate" % (en+1,len(keys)))
	st = snp_id - rng
	ed = snp_id + rng
	if st < pass_n:
		st = pass_n
	if ed > len(vcf):
		ed = len(vcf)
	print("[PROCESS] VCF Slice mode activated")	
	print("[PROCESS] %s SNPs are target to assess genotype" % (ed-st+1))
else:
	print("[PROCESS] VCF Slice mode disabled, whole variants would be output.")

out = open(prefix,"w")
if mode_slice == True:
	for i in vcf[st:ed]:
		i = i.split()
		idx = i[2]
		header.append(idx)
		for en,sp in enumerate(i[9:]):
			sp = sp.split(":")[0].split("/")
			try:
				sp = int(sp[0]) + int(sp[1])
			except ValueError:
				sp = "NA"
			geno[keys[en]].append(sp)

	print(*header,sep="\t",file=out)
	for key in keys:
		print(key,end="\t",file=out)
		print(*geno[key],sep="\t",end="\n",file=out)
else: ## mode_slice == False:
	print(*header,sep="\t", file=out)
	keys = sorted(list(geno.keys()))
	for k in keys:
		print(k, end="\t", file=out)
		print(*geno[k], sep="\t", file=out)
		
out.close()	
		
