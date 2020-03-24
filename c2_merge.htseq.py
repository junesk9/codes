#!/usr/bin/env python3


"""
Merging multiple HTseq-count output files to single file,
the output will be STDOUT.


AUG 23rd 2016 Junesk9


[Change log]
Feb 08th 2017	Add new options "-d" and "-r", for aligning samples by numeric indicator at middle
Dec 26th 2016	Add new option "-o" for alternative output file name
		Add progress dialogs

"""

import sys
import glob

out = {"-o":"merged.HTseq.txt", "-d":"_", "-r":2, "sort":False}
args = sys.argv[1:]
if "-o" in args:
	out_idx = args.index("-o")+1
	out["-o"] = args[out_idx]
if "-sort" in args:
	out["sort"] = True
if "-d" in args:
	del_idx = args.index("-d")+1
	out["-d"] = args[del_idx]
if "-r" in args:
	rank_idx = args.index("-r")+1
	try:
		out["-r"] = int(args[rank_idx])
	except:
		print("[ERROR] The rank option ("-r") accept an integer only")
		sys.exit(2)
else: pass

sort = out["sort"]
dlm = out["-d"] ## delimeter
rnk = int(out["-r"])-1 ## ranker

out = out["-o"] if out["-o"].endswith("HTseq.txt") else out["-o"] + ".HTseq.txt"



fi_list = glob.glob("*txt")
fi_list = sorted(fi_list)
if sort == True:
	print("[PROGRESS] Sorting option is given, delimeter: %s" % dlm)
	fi_list = sorted(fi_list, key=lambda x: int(x.split(dlm)[rnk]))


print("[PROGRESS] %s COUNT files detected" % len(fi_list))
print("[PROGRESS] The output file will be %s" % out)



dic = {}
header = ["ID"]
#print(fi_list)
#sys.exit()

for fi in fi_list:
	head = fi.split(".")[0]
	header.append(head)
	for i in open(fi):
		i = i.strip().split()
		#print(i)
		#sys.exit()
		idx = i[0]
		val = i[1]
		try:
			dic[idx] += [int(val)]
		except KeyError:
			dic[idx] = [int(val)]

keys = sorted(list(dic.keys()))

#print(out)
#print(header)

out = open(out,"w")
print(*header,sep="\t",file=out)
for key in keys:
	if sum(dic[key]) != 0 and not key.startswith("_"):
		print(key,end="\t",file=out)
		print(*dic[key],sep="\t",file=out)
out.close()
print("[PROGRESS] All progress done !")
print()
