#!/usr/bin/env python3 

"""
Assign peaks from the BedGraph files of sRNA mapping

Continous read depth (more than 'cutoff') will be merged into single peak
the assigned peak number will be tagged at fifth additional column of Bed
[chr][st][ed][depth][peak#]

[USAGE]

python3 s5_bedPeakMaker.py (--count)

two Output files made for each BedGraph file:
1. A Bedfile contains only peak assigned lines
2. Peak info list files - Adaptable as input for 'm5_find.DMR.prox.gene.py'
[peak][chr][st][ed][width][depth] 

*Count-only mode added (For time-saving; with already made [1]outfile)
if --count added in the command, this will only run count modes with bedfiles in the current folder.
peak length is adjustable, but peak depth.
the output will be made in the current folder

15. 02. 12. Toshimori

Change logs--
15.03.28 Count-only mode added. 
15.03.28 Add new clumn for 2nd output file [depth], stand for total nulceotides mapped to single peak

"""


import glob, os, time, sys

## Swiching mode count-only or not 
args = sys.argv[1:]
if "--count" in args:
	count = "yes"
	print("### This run with the count-only mode")
else: count = "no"

#---------------------------Arguments
cutoff = 10 # read depth restriction for forming a peak
width = 50 # width restriction for peak width
inputls = glob.glob("*bedgraph")

### Make output folder
stime=time.time()
if count == "no":
	outfolder = "peak_"+str(int(stime))
	os.mkdir(outfolder)
else: 
	print("This run with both making peak.bed and count modes")


## Make peak-annotated bed file
for fi in inputls:
	if count == "no":
		prefix = outfolder+"/"+(".").join(fi.split("/")[-1].split(".")[:2])+".more"+str(cutoff)+".peak.bedgraph"
		prefix2 = outfolder+"/"+(".").join(fi.split("/")[-1].split(".")[:2])+".more"+str(cutoff)+"_"+str(width)+"nt"+".peak.txt"
		out = open(prefix,"w")
		k,end = 0,0
		for j in open(fi):
			j = j.strip()
			height = float(j.split("\t")[3]) #scientific notation (e.g. 1.0+10E6) is not recognized by "integer" but "float"
			st = int(j.split("\t")[1])
			ed = int(j.split("\t")[2])
			#print(st,ed,end)
			if height >= cutoff and st != end:
				k += 1
				peak = "peak"+str(k)
				end = int(j.split("\t")[2])
				print(j,peak,sep="\t", file=out)
			if height >= cutoff and st == end:
				peak = "peak"+str(k)
				end = int(j.split("\t")[2])
				print(j,peak,sep="\t", file=out)
			else: pass
		out.close()
		ctime = time.time()
		print("Make peak.bedgraph Done! %s, takes %s secs" % (fi,int(ctime-stime)))
	else: pass

	## Different file naming for count-only mode
	if count == "yes":
		prefix2 = ".".join(fi.split("/")[-1].split(".")[:2])+".more"+str(width)+"nt"+".peak.txt"
		fi2 = fi
	else: 
		fi2 = prefix
	#end of if
	
	##--------------------------------body for make list of sRNA peak and count read-depth
	out2 = open(prefix2, "w")
	print("#peak","chr","peak_st","peak_ed","peak_len","depth",sep="\t",file=out2)
	dic = {}
	for k in open(fi2):
		k = k.strip().split("\t")
		peak = k[-1]
		ch = k[0]
		st = int(k[1])
		ed = int(k[2])
		depth = float(k[3])
		try:
			if st < dic[peak][1]:
				dic[peak][1] = st
			if ed > dic[peak][2]:
				dic[peak][2] = ed
			else: pass
			dic[peak][3].append(depth)
		except:
			dic[peak] = [ch,st,ed,[depth]]
		
	keys = sorted([i for i in dic.keys()])
	for x in keys:
		peakLen = int(dic[x][2]) - int(dic[x][1])
		dic[x][3] = str(int(sum(dic[x][3]))) ## the integer sum of the list of float numbers
		if peakLen >= width:
			dic[x].insert(0, x)
			dic[x].insert(-1,str(peakLen))
			for i in dic[x]:
				print(i,end="\t",file=out2)
			print(file=out2)
		else: pass
	out2.close()
	ctime = time.time()
	print("Make peak.list Done! %s, takes %s secs" % (fi,int(ctime-stime)))
dic, keys = [],[] ## Make the memory blank
