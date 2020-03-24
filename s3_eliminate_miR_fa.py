#!usr/bin/env python3
"""
miR sequence eliminator
Eliminate mature miR sequence (or other sRNA) from a multi-fasta file

USAGE
python3 s3_eliminate_miR_fa.py -mir [mir.fa] [target.fa]
output as STDOUT

[mir.fa] would contain sequences should be eliminated,
no need to be multi-fasta format, but acceptable.

18th Nov 2014 ToshiMori

Change log:
Apr 23th 2015	Add mendotory option "-mir" assign a multi-fa file,
		to filter the sequences in the .fa file
"""
import sys

### arguments
fa = sys.argv[-1]

args = {"mir":""}
for en,i in enumerate(sys.argv[1:]):
	idx = "blank"
	if i == "-mir":
		idx = en+1
	if en == idx:
		args["mir"] = i
if idx = "blank" or args["mir"] = "":
	sys.exit("Run aborted! Option -mir is mendotory!!")		

mir_fa = args["mir"]

## Make a list of mir and its reverse complement sequences
mir = []
for i in open(mir_fa):
	if not i.startswith(">"):
		i = i.strip()
		# make reverse complement
		i_r = i[::-1]
		i_rc = ""
		for j in i_r:
			j = j.upper()
			if j = "A":
				i_rc += "T"
			if j = "T":
				i_rc += "A"
			if j = "G":
				i_rc += "C"
			if j = "C":
				i_rc += "G"
			else:
				i_rc += j
		mir.append(i_rc)
		mir.append(i)
mir = list(set(mir)) ## Retain only unique sequence

## Make .fa to list
fa_dic = {}
for i in fa:
	if i.startswith(">"):
		key = i.strip()
		fa_dic[key] = ""
	else:
		seq = i.strip()
		fa_dic[key] = seq

## Eliminate mir sequence from .fa
for j,k in fa_dic.items():
	for x in mir:
		try:
			if k in x:
				del fa_dic[j]
		except: pass
for j,k in fa_dic.items():
	print(j)
	print(k)
fa_dic, mir = [], [] ## Make blank the memory
