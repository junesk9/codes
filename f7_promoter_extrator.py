#!/usr/bin/env python3

"""
Promoter extractor from Genome with GFF3 files
Only promoters with the length given output


[USAGE]
python3 promoter_extractor.py -genome [] -gff []
-genome: multi-fasta formatted genome file
-gff:	 GFF3 formatted genome annt file

[OPTIONS]

-len:	Restrict the length of the promoter, [default: -1000]
-out:	Determine the output file name
-gene:	Restrict genes to be collected, 
	a flat txt file containing gene names [default: whole genes in GFF]
-negl:  Neglect the input gene list, the output contains results from all genes but from input list
        [default: False]

==CHANGE LOG
18JAN 2017	Minor upgrades, add option "-negl"
28OCT 2015	First build
"""

import sys

################### Arguments

opt = {"-genome":"","-gene":"","-len":1000,"-out":"","-gff":"", "-negl":False}
args = sys.argv[1:]

try:
	genome_idx = args.index("-genome")+1
	opt["-genome"] = args[genome_idx]
	gff_idx = args.index("-gff")+1
	opt["-gff"] = args[gff_idx]
except:
	print()
	print("======= Promoter extractor =========")
	print("[ERROR]Argument for [-genome] and/or [-gff] is missing!")
	print(__doc__)
	sys.exit(2)

if "-gene" in args:
	gene_idx = args.index("-gene")+1
	opt["-gene"] = args[gene_idx]
if "-len" in args:
	len_idx = args.index("-len")+1
	try:
		opt["-len"] = int(args[len_idx])
	except:
		print("The [-len] option accepts integer only!!")
		print(__doc__)
		sys.exit(2)
if "-negl" in args:
	opt["-negl"] = True

genome, gene, length, gff, out, negl = opt["-genome"],opt["-gene"],opt["-len"],opt["-gff"],opt["-out"],opt["-negl"]


if out == "":
	if gene == "":
		out = genome.split(".")[0]+"."+str(length)+"nt.promoter.fa"
	else:
		out = gene.split(".")[0]+"."+str(length)+"nt.promoter.fa"
else: pass

print("""

	########Promoter extractor########
	
	genome: %s 
	gff   : %s
	genes : %s (default: blank(all))
	length: %s
	out   : %s
	negl  : %s
	""" % (genome, gff, gene, length, out, negl))


########################## Methods
def rev_com(nt):
	tmp = ""
	for i in nt[::-1]:
		i = i.upper()
		if i == "A":
			tmp += "T"
		elif i == "T":
			tmp += "A"
		elif i == "C":
			tmp += "G"
		elif i == "G":
			tmp += "C"
		else:
			tmp += i
	"""
	if len(nt) != len(tmp):
		print(len(nt),len(tmp))
		print(nt)
		print(tmp)
	"""
	return tmp


#################### Load genome

tmp = {}
for i in open(genome):
	i = i.strip().split()[0]
	if i.startswith(">"):
		name = i.split(".")[0][1:]
		tmp[name] = []
	else:
		tmp[name].append(i)

for i in tmp.keys():
	tmp[i] = "".join(tmp[i])

genome = tmp


print("[STEP1] Genome fasta is loaded")
################### Load GFF


tmp = {} # gene_id : ["chr","ori","st","ed"]
for i in open(gff):
	i = i.strip().split("\t")
	#print(i)
	#sys.exit()
	if i[0].startswith("#"):
		pass
	elif len(i) > 2 and i[2] == "gene":
		ch = i[0]
		ori = i[6]
		st = int(i[3])
		ed = int(i[4])
		g_id = i[-1].split(";")[0].split(":")[-1].upper() ## prevent confusions in upper/lower cases
		tmp[g_id] = [ch,ori,st,ed]
gff = tmp


print("[STEP2] %s of Gene info is loaded" % len(gff.keys()))
################### if gene_list was provided

print(gene)

if gene == "":
	g_list = list(gff.keys())
else:
	g_list = [i.strip().split(".")[0] for i in open(gene)]


for i in g_list:
	if i.startswith(">"):
		naked = i[1:]
		g_list.remove(i)
		g_list.append(naked)
#end of for
g_list = sorted(list(set(g_list))) ## Remove duplications & sorting by names



### Process "negl" option
tmp = list(gff.keys())
x, y = 0, 0
if negl == True:
	for g in g_list:
		y += 0
		gi = g.upper()
		if gi in tmp:
			x += 1
			tmp.remove(gi)	
		else: pass
	#end of for
	print("[STEP3] the Negl option activated, %s attempted, %s removed from %s total genes" % (y,x,len(gff.keys())))
	g_list = tmp
else: pass
##### remove ambigous (not included in GFF) genes from gene list
x, y = 0, 0
for g in g_list:
	gi = g.upper()
	if gi in gff.keys():
		pass
	else:
		g_list.remove(g)
		x += 1
		print("[CAUTION] Gene ID %s was not included in the genome" % g)
	y += 1
print("[STEP3] %s gene IDs attempted, %s IDs retracted" % (y,x))

#################### Extract promoter seq

count = 0
out = open(out,"w")
for i in g_list:
	idx = i.upper() ## Some GFFs have gene names of upper cases.
	ch = gff[idx][0]
	ori = gff[idx][1]
	st = gff[idx][2]-1
	ed = gff[idx][3]
	if ori == "+":
		rng = [st-length,st]
	elif ori == "-":
		rng = [ed,ed+length]
	#end of if
	if rng[0] < 0: ## at the left tip of a chromosome
		rng[0] = 0
	#end of if
	nt = genome[ch][rng[0]:rng[1]]
	if ori == "-": 
		#print(i,len(nt))
		nt = rev_com(nt)
		#print(i,len(nt))
	"""
	if i == "Traes_6BS_FB06D4612":
		print(ch,ori,st,ed,rng)
		print(nt)
		sys.exit(2)
	"""
	#print(i,len(nt),ori,rng,len(genome[ch]))
	if len(nt) == abs(length):
		count += 1
		ori = "fwd" if ori == "+" else "rev"
		logo = ">" + i + ".prom" + str(length) +"." + ori
		print(logo,file=out)
		print(nt,file=out)
print("[STEP4] %s/%s of promoter regions are printed " % (count,len(gff.keys())))
out.close()
