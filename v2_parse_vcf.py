#!usr/bin/env python3

"""
VCF file parser with genetic information
Manipulate single VCF file to simplified snp.txt
[gene][site][chr][gene_st][gene_ed][ref][allele][allele_count][allele_count][frequency][snp_qual]
# allele_count adopt from [DP4] value in VCF
# SNP_quality  adopt from SNP quality score at VCF file (column 6)
# Genetic information adopt from any GFF3 file


May, 28th, 2015 Toshimori
"""

import glob, sys


fi_ls = glob.glob("../../2_vcf/*vcf")
info = open("Radish_V1.51.gff")
info_dic = {}
for x in info.readlines():
	x = x.strip().split()
	if x[2] == "gene":
		gene = x[-1].split("=")[-1]+".1"
		ch = x[0]
		st = x[3]
		ed = x[4]
		ori = x[6]
		info_dic[gene] = [ch,st,ed,ori]


for fi in fi_ls:
	prefix = ".".join(fi.split("//")[-1].split(".")[:-1])+".snp.txt"
	out = open(prefix, "w")
	print("gene","site","chr","gene_st","gene_ed","ref","allele","allele_count","ref_count","freq","snp_qual",sep="\t",file=out)
	for i in open(fi).readlines():
		if not i.startswith("#"):
			i = i.strip().split("\t")
			gene = i[0]
			site = i[1]
			ref = i[3]
			snp = i[4]
			qual = i[5]
			annt = info_dic[gene]
			if "DP4" in i[7]:
				info = i[7].split(";")
				for j in info:
					if j.startswith("DP4"):
						j = j.split("=")[-1].split(",")
						k_ref = int(j[0])+int(j[1])
						k_snp = int(j[2])+int(j[3])
						k_rate = round(k_snp/(k_snp+k_ref),2)
			else: 
				k_snp = "na"
				k_ref = "na"
				k_rate = "na"
			print(gene,site,annt[0],annt[1],annt[2],ref,snp,k_snp,k_ref,k_rate,qual,sep="\t",file=out)

	out.close()
