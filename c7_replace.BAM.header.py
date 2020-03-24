#!/usr/bin/env python3

"""
Automated BAM header changer 


2018.11.19 Junesk9


"""

import sys
from subprocess import call
import glob


#############Args

RG = {"ID":"group1",
	"SM":"null",
	"PL":"illumina",
	"LB":"lib1",
	"PU":"unit1"}



############## Running body
bams = glob.glob("*bam")
tot = len(bams)
print("[PROGRESS] Detected %s BAM files" % tot) 

x = 0
for bam in bams:
	sm = bam.split(".")[0]
	RG["SM"] = sm
	h_rg = "@RG	ID:%s	SM:%s	PL:%s	LB:%s	PU:%s" % (RG["ID"],RG["SM"],RG["PL"],RG["LB"],RG["PU"])
	#print(h_rg)

	command = "samtools view -H %s > old.header" % bam
	call(command, shell=True)
	
	h_out = open("new.header","w")
	for h in open("old.header"):
		h = h.strip()
		if h.startswith("@RG"):
			print(h_rg,file=h_out)
		else:
			print(h, file=h_out)
	h_out.close()

	newbam = bam[:-3]+"rehead.bam"
	command = "samtools reheader new.header %s > %s" % (bam, newbam)
	cmd2 = "rm *header"
	cmd3 = "rm %s" % bam
	call(command, shell=True)
	call(cmd2, shell=True)
	call(cmd3, shell=True)

	x += 1
	print("[PROGRESS] Progress %s Done [%s/%s]" % (newbam, x, tot))
