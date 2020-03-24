#!/usr/bin/env python3

"""
Generate HOMER peak file from AGI list (both transcript & gene_ID are applicapable)

Since in present (May2017) HOMER cannot process the AGI list (findMotif.pl)
so temperaly convert the AGI code to HOMER peak to apply (findMotifGenome.pl).

>generate.HOMER-peak.py [input] -> output.txt

>sudo su ## dont know why but need to be SU
>source /home/alrcpmb1/.profile

>findMotifGenome.pl output.txt arabidopsis homer-out-folder -len 6,8 -p 20


16th MAY 2017 Junesk9


"""

import sys 


opt = {"inp":"", "gtf":"araport11.gtf"}
opt["inp"] = sys.argv[-1]


################ Argument
prom_len = 600 # promoter length from CDS start site

################ Parsing gene list
inp = opt["inp"]
inp = [i.strip() for i in open(inp)]
inp = [i[:9] if "." in i else i for i in inp]
##print(inp[:5])



############## Parsing GTF
### peak line [ID, CHR, ST, ED, STRADN (+/- or 0/1)]
gtf = opt["gtf"]
gtf = [i.strip().split("\t") for i in open(gtf)]

tmp = {}
tmp2 = {}
for g in gtf:
	"""
	if g[2] == "start_codon":
		gid = g[-1].split()
		gdx = gid.index("gene_id")
		#gid = gid[gdx+1][1:][:-1]
		tdx = gid.index("transcript_id")
		tid = gid[tdx+1][1:][:-2]
		ch = g[0][-1]
		st = g[3]
		ed = g[4]
		ori = g[6]
		if ori == "+":
			tp = [tid,ch,st,ed,ori]
		else:
			tp = [tid,ch,ed,st,ori]
		tmp[tid] = tp
	"""
	if g[2] == "exon":
		gid = g[-1].split()
		gdx = gid.index("gene_id")
		gid = gid[gdx+1][1:][:-1]
		ch = g[0][-1]
		rng = [int(g[3]), int(g[4])]
		ori = g[6]
		try:
			tmp[gid] += rng
		except:
			tmp[gid] = rng
			tmp2[gid] = [ch,ori]

tmp3 = {}
for t in tmp.keys():
	st = min(tmp[t])
	ed = max(tmp[t])
	ch = tmp2[t][0]
	ori = tmp2[t][1]
	line = [t,ch,st-prom_len,st-1,ori]
	tmp3[t] = line

		
for i in inp:
	print(*tmp3[i],sep="\t")
		
