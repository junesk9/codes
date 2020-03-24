#!/usr/bin/env python3


"""
UPDATED NOV 19th 2018 Junesk9

Now available with HTSeq-count outputs
Now also process Kallisto outputs.
The usage is same to that for eXpress,
1. Generate list of outputs (e.g. abundance.tsv) as a single txt
2. run the command with the txt generated above
3. the code detect the outputs from whether eXpress, Kallisto, or HTSeq-count, and merge the data to tsv files.
4. the lines with value sum of ZERO go discarded.

[USAGE]
c6_merge.express.results.py [output's list]


Jun 07th 2017 Junesk9

Command to merging eXpress results from RNA-seq analysis
outputs 2 files contains "tot_cnt" or "fpkm" as columns of samples
the lines with value sum of 0 goes discarded


"""

import sys

ls = sys.argv[-1]
if len(sys.argv) == 1:
	ls = "result.txt"
ls = sorted([i.strip() for i in open(ls)])


if  ls[0].endswith("results.xprs"):
	mode = "EXPRESS"
elif ls[0].endswith("abundance.tsv"):
	mode = "KALLISTO"
elif open(ls[0]).readlines()[-1].startswith("__alignment_not_unique"):
	mode = "HTSEQ"
else:
	print("[ERROR] Unknown output types")
	print(ls[0])
	sys.exit()
print("[PROGRESS] Detected %s outputs" % mode)

if mode == "EXPRESS":
	cnt, fpkm = {}, {}
	header= ["transcript"]


	for l in ls:
		h = l.split("/")[1]
		header.append(h)
		for i in open(l.strip()):
			if i.startswith("bundle_id"):
				pass
			else:
				i = i.strip().split()
				gid = i[1]
				c = int(i[4])  ## Total Count
				f = round(float(i[10]),3) ## FPKM
				try:
					cnt[gid].append(c)
					fpkm[gid].append(f)
				except KeyError:
					cnt[gid] = [c]
					fpkm[gid] = [f]

	outc = open("eXpress.count.merged.txt","w")
	outf = open("eXpress.fpkm.merged.txt","w")
	print(*header,sep="\t",file=outc)
	print(*header,sep="\t",file=outf)
	keys = sorted(list(cnt.keys()))
	for k in keys:
		if k.endswith(".1") or k.endswith("GFP"):
			if sum(cnt[k]) != 0:
				print(k,end="\t",file=outc)
				print(*cnt[k],sep="\t",file=outc)
			if sum(fpkm[k]) != 0:
				print(k, end="\t",file=outf)
				print(*fpkm[k],sep="\t",file=outf)
	outc.close()
	outf.close()

elif mode == "KALLISTO":
	cnt, tpm = {}, {}
	header = ["Transcript"]

	for l in ls:
		h = l.split("/")[1]
		header.append(h)
		for i in open(l):
			if i.startswith("target_id"):
				pass
			else:
				i = i.strip().split()
				gid = i[0]
				c = float(i[3]) ## Est_counts
				f = float(i[4]) ## tpm (transcript per million)
				try:
					cnt[gid].append(c)
					tpm[gid].append(f)
				except KeyError:
					cnt[gid] = [c]
					tpm[gid] = [f]

	#print(list(cnt.keys())[:10])
	#sys.exit()


	outc = open("Kallisto.count.merged.tsv","w")
	outf = open("Kallisto.tpm.merged.tsv","w")
	print(*header, sep="\t", file=outc)
	print(*header, sep="\t", file=outf)
	keys = sorted(list(cnt.keys()))
	for k in keys:
		if sum(cnt[k]) != 0:
			print(k,end="\t",file=outc)
			print(*cnt[k],sep="\t",file=outc)
		if sum(tpm[k]) != 0:
			print(k, end="\t", file=outf)
			print(*tpm[k], sep="\t", file=outf)
	outc.close()
	outf.close()
elif mode == "HTSEQ":
	cnt = {}
	header = ["Gene_id"]

	for l in ls:
		h = l.split("/")[1].split(".")[0]
		header.append(h)
		for i in open(l):
			if i.startswith("__"):
				pass
			else:
				i = i.strip().split()
				gid = i[0]
				c = int(i[1])
				try:
					cnt[gid].append(c)
				except KeyError:
					cnt[gid] = [c]

	out = open("HTSEQ.count.merged.tsv","w")
	print(*header, sep="\t", file=out)
	keys = sorted(list(cnt.keys()))
	for k in keys:
		if sum(cnt[k]) != 0:
			print(k, end="\t", file=out)
			print(*cnt[k], sep="\t", file=out)
	out.close()
