#!/usr/bin/env python3


"""
BAM to FASTQ 

Equips both decomprassing GZIP and Loading Binary BAM 
(also possible with plain TXT FASTQ and SAM; automatically determined)

Two modes are prepared:
I. Extract FASTQ from the BAM/SAM file
 -bam option required
 the fastq information extracted from the input BAM/SAM file 
 & output to the [prefix].fastq plain text file.

II. Slice input FASTQ owing to input BAM
  Both inputs are required;
  -bam : input bam or sam file
  -fq  : input fastq file (both gZIP and plain fastq)

  The output STDOUT FASTQ format


--- log ---
2017.08.01 De-bugging; Add the extract mode
2017.07.31 First build

"""

import sys
import mimetypes ## To detect filetype
import subprocess ## Popen() well, not easy to understand,
import gzip, bz2


####################### Arguments
args = {"bam":"","fq":""}

arg = sys.argv[1:]
#print(arg)
if "-bam" in arg:
  b_index = arg.index("-bam")+1
  args["bam"] = arg[b_index]
if "-fq" in arg:
  a_idx = arg.index("-fq") + 1
  args["fq"] = arg[a_idx]

#INPUT validation
if args["bam"] == "":
  print("[ERROR] -bam is mendatory !! ")
  print()
  sys.exit()
#end of if


bam,fq = args["bam"],args["fq"]

# Determine mode whether referring the external FASTQ file or not
extr_mode = True if fq == "" else False

#Validate the input file as SAM or BAM
mime = mimetypes.guess_type(bam)[0]
bam_mode = True if mime == None else False


######################## BAM/Manipulation

def run_and_capture(cmd):
  """
  param cmd: str command
  rtype: str
  retun: stdout
  """

  proc = subprocess.Popen(cmd, shell=True,  stdout=subprocess.PIPE)
  buf = []
  z = 0  

  for line in proc.stdout.readlines():
    z += 1
    #print(line[0].strip())
    #sys.exit()
    #print(line)
    line = line.decode("utf-8").split() ## Convert Bytes to Strings
    #print(z, line)
    buf.append(line)
    #sys.stdout.write(line)
    #print([proc.poll()])

    if not line and proc.poll() is not None:
      break

  proc = [] ## Flush the memory
  return buf

########################### Slice FASTQ referring BAM/SAM file
if extr_mode == False:
  if bam_mode == True:
    cmd = "samtools view " + bam 
    sam = run_and_capture(cmd)
    #print(sam[0])
    """ 
    except:
      print("[ERROR] SAMTOOLS running error !")
    """
  else: #bam_mode == False:
    sam = [i.strip().split() for i in open(bam)]


  # Extract fastq ids
  tmp = []
  for i in sam:
    ids = i[0]
    ids = ids[1:] if ids[0] == "@" else ids
    tmp.append(ids)
  sam = sorted(list(set(tmp)))
  #print(len(sam))
  #sys.exit()

  ####################### FASTQ manipulation

  t = mimetypes.guess_type(fq)[1]
  gz = True if t == 'gzip' else False

  if gz == True:
    fq = gzip.open(fq,mode="rt")

  fq_lib = {}
  for en,i in enumerate(fq):
    if en % 4 == 0:
      name = i.split()[0][1:] ## Discard the first "@" as the seq id indicator
      fq_lib[name] = []
    else:
      i = i.strip()
      fq_lib[name].append(i)
  fq = fq_lib

  ###################### Search and Generate output
  for n in sam:
    try:
      seq = fq[n]
    except KeyError:
      print("[ERROR] Cannot find the seq name of %s" % n)
      print()
      sys.exit()
    name = "@"+n
    print(name)
    print(*seq,sep="\n")

############################## Extract FASTQ from BAM/SAM file
elif extr_mode == True:
  print("[PROGRESS] FASTQ extracted from %s" % bam)
  out = ".".join(bam.split(".")[:-1])+".fastq"
  print("[PROGRESS] FASTQ saved into %s" % out)
  if bam_mode == True:
    cmd = "samtools view " + bam
    sam = run_and_capture(cmd)
  else:
    sam = [i.strip().split() for i in open(bam)]
    
  out = open(out, "w")
  for s in sam:
    #print(s)
    if s[0].startswith("@"):
      pass
    else:
      f1 = "@" + s[0]
      f2 = s[9]
      f3 = "+"
      f4 = s[10]
      print(f1,f2,f3,f4,sep="\n",file=out)
  out.close()




