#!/usr/bin/env python3

"""

USAGE:
python3 shirink_vcf_RNA-seq.py [input_VCF]
 

"""

import sys

args = sys.argv[1:]

in_vcf = args[0]
prefix = ".".join(in_vcf.split(".")[:-1])
ref_csv = prefix+".annt.csv"
out = prefix+".shirinked.csv"

opt = {"-csv":ref_csv,"-cds_only":True,"-out":out}

if "-csv" in args:
	idx_csv = args.index("-csv")
	opt["-csv"] = args[idx_csv+1]
if "-out" in args:
	idx_out = args.index("-out")

ref_csv, cds_only, out = opt["-csv"], opt["-cds_only"], opt["-out"]

print("""
	######## CSV selector with anntated CSV data #####
	
	This program runs with following options:
	-input_vcf    : %s 
	-ref_csv      : %s
	-out_file     : %s
	-cds_only     : %s

""" % (in_vcf, ref_csv, out, cds_only))


################################## Load SNP annotation data (CSV format)
try:
	ref_csv = [i.strip().split(",") for i in open(ref_csv)]
except:
	print("[ERROR] cannot locate the csv annotation file: %s" % ref_csv)
	print(__doc__)
	sys.exit()

ref = []
for i in ref_csv:
	info = i[2]
	snp = i[1]
	if "Effective" in info or "Silence" in info:
		ref.append(snp)
	elif cds_only == False:
		if "intronic" in info:
			ref.append(snp)

print("Loading SNP anntation data Done, Totally %s SNPs, filtered %s SNPs\t" % (len(ref_csv),len(ref))
################################# Load VCF and parse it refering the CSV annotation data

try:
	in_vcf = open(in_vcf)
except:
	print("[ERROR] cannot locate the input vcf file: %s" % in_vcf)
	print(__doc__)
	sys.exit()


header, collected = [], []
x,y,z = 0,1,10000
for i in in_vcf:
	i = i.strip()
	if i.startswith("#"):
		header.append(i)
	else:
		snp = i.split("\t")[2]
		if snp in ref:
			collected.append(i)
	x += 1
	if x >= z:
		print("Processing %s lines done!" % (y*z))
		x = 0
		y += 1

############################### Make output file

print("Make output CSV file: %s, containing %s snp info" % (out, len(collected)))

out = open(out,"w")
print(*header,sep="\n",file=out)
print(*collected,sep="\n",file=out)
out.close()
