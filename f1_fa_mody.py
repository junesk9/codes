#!/usr/bin/env python3
"""

### Change log
15.01.29	Modify the FQ2FA module
14.12.12	Equip faster mulit-fa loading tool
"""
import sys
import os
import time
import gzip, bz2
   

args = sys.argv[1:]

##############################################################
###### Functions

def fq2fa(sFastafile_open):
    """
    convert FastQ to FastA
    """

    for n, sLine in enumerate(sFastafile_open):
        if n % 4 == 0 and sLine.startswith("@"):  # The 1st line of FQ file must start with "@".
            print(">"+sLine[1:].strip().split(" ")[0])
        #end of if
        if n % 4 == 1:
            if sLine.startswith("N"):
                sLine = sLine[1:]
            print(sLine.strip())
        #end of if
    #end of for

def Make_DataSet_fastaData(sFastafile_open):
    """
    Convert multi-fasta file to Dictionary db.
    """

    dFastaDataSet={}
    key,seq = [],[]
    c = -1
    for sLine in open(args[-1]):
        if sLine.startswith('>'):
            fasta_key = sLine.strip().split(' ')[0]
            fasta_key = fasta_key[1:]
            key.append(fasta_key)
            seq.append([])
            c += 1
        else:
            ss = sLine.strip()
            seq[c].append(ss)


    for e,i in enumerate(seq):
        i = "".join(i)
        seq[e] = i

    for e,j in enumerate(key):
        dFastaDataSet[j] = seq[e]
    return dFastaDataSet


def Write_split_fasta(dFastaDataSet):
    """
    make single fasta files and store at ./[filename]_out
    """

    sPathName = args[1].split(".")[0]   # Collect file name as output folder name
    if os.path.exists(sPathName+"_out"):  # Verbose for making output folder
        print("== path ["+sPathName+"_out] already exist ==")
    else:
        print("== make path [" + sPathName+"_out] ==")
        os.mkdir(sPathName+"_out")
    #end of if
    os.chdir("./"+sPathName+"_out")  # Making output folder
    print("== outputs stored at ["+sPathName+"_out] ==\n")
    
    for sKey, sValue in dFastaDataSet.items():  # Collapse the fasta name as output file name, collect first word only
        sKey_ls = sKey.replace(" ", "|").replace("_", "|").split("|")
        sName = sKey_ls[0]
        sFilename = sName+'.fa'

        print(sFilename)
        with open(sFilename, mode='w') as out:
            out.write('>'+sKey+'\n'+sValue)
        #end of with
    #end of for

def Calculate_fasta_stat(dFastaDataSet):
    """
    Calculate sequence number, total basepairs,
    average basepairs, N50-N90, and ATGCN contents
    """
    
    nLine_count = len(dFastaDataSet.values())  # Calculte length stats
    # nBp_count = 0 # for calculate tot. length, but not used; subs to len(sSeq)
    nLen_ls = []  # for N50 calculation
    sSeq_ls = []  # for calculate length & ATGCN contetents
    for i in dFastaDataSet.values():
        ##nBp_count += len(i)
        nLen_ls.append(len(i))
        sSeq_ls.append(i)
    #end of for
    sSeq = "".join(sSeq_ls)
    nLen_ls = sorted(nLen_ls, reverse=True)
    print("\n\nNumber of seq:",nLine_count)
    print("Total bps:", len(sSeq))
    print("Average bps:", round(len(sSeq)/nLine_count,2))
    print("Maximum bps:", nLen_ls[0])
    print("Minimum bps:", nLen_ls[-1])
           

    a = 1                                     # Calculate N50 and further
    b = 0
    nSum = 0
    for j in nLen_ls:
        b += 1
        nSum += j
        if nSum >= len(sSeq)*0.5 and a == 1:
            print("N50: %s\t# %s" % (j, b))
            a += 1
        if nSum >= len(sSeq)*0.6 and a==2:
            print("N60: ", j)
            a += 1
        if nSum >= len(sSeq)*0.7 and a==3:
            print("N70: ", j,"\t# ",b)
            a += 1
        if nSum >= len(sSeq)*0.8 and a==4:
            print("N80: ", j)
            a += 1
        if nSum >= len(sSeq)*0.9 and a==5:
            print("N90: ", j)
            a += 1
        
        #end of if
    #end of for

    '''[[old ver]]
    ATGC_dic = {"A":0,"T":0,"C":0,"G":0,"N":0} ##Calculate ATGCN contents
    try:
        for k in sSeq:
            ATGC_dic[k] += 1
        #end of for
    except Exception as e:  ## If error raised by other character found from sSeq, raise error msg & continue process
        print("%r" % e)
    #end of try

    print("A: ",ATGC_dic["A"]," (",round(ATGC_dic["A"]/len(sSeq)*100,2),"%)",sep="", end="\t")
    print("T: ",ATGC_dic["T"]," (",round(ATGC_dic["T"]/len(sSeq)*100,2),"%)",sep="", end="\t")
    print("G: ",ATGC_dic["G"]," (",round(ATGC_dic["G"]/len(sSeq)*100,2),"%)",sep="", end="\t")
    print("C: ",ATGC_dic["C"]," (",round(ATGC_dic["C"]/len(sSeq)*100,2),"%)",sep="", end="\t")
    print("N: ",ATGC_dic["N"]," (",round(ATGC_dic["N"]/len(sSeq)*100,2),"%)‰",sep="", end="\t")
    '''

    Str_dict = {}  # Calculate ATGCN contents
    for k in sSeq:
        try:
            Str_dict[k.upper()] += 1  # If k as key present in dict, count-up
        except:
            Str_dict[k.upper()] = 1  # or add k as a new key
        #end of try
    #end of for

    for l, m in sorted(Str_dict.items()):   # load dict and print the data
        perc = round(m/len(sSeq)*100, 2)  # calculate percentage
        print('{0}: {1} ({2}%)'.format(l, m, perc), end="\t")  # use '.format'
    #end of for

def Split_single_to_multi(dFastaDataSet):
    """
    For both multi-fasta and single-fasta file,
    split fasta sequence to multi-lines with 60 line length
    """

    nLine_len = 60  # Change '60' for the initiative line length

    for sKey, sValue in dFastaDataSet.items():
        print(">"+sKey)
        for i in range(0, len(sValue), nLine_len):
            print(sValue[i:i+nLine_len])
        #end of for
    #end of for
 

##############################################################
###### Main            


def main(args):
    cTime = time.time()

    if len(args) == 2 and args[0] == 'FQ2FA':
        sFastafile_open = open(args[1], 'rU')
        sys.exit(fq2fa(sFastafile_open))

    elif len(args) == 2 and args[0] == 'MF2SF':
        sFastafile_open = open(args[1], 'rU')
        dFastaDataSet = Make_DataSet_fastaData(sFastafile_open)
        Write_split_fasta(dFastaDataSet)
        sFastafile_open.close()
        print("\nDone! it takes", round(time.time()-cTime, 4), "sec")

    elif len(args) == 2 and args[0] == 'STAT':
        sFastafile_open = open(args[1], 'rU')
        dFastaDataSet = Make_DataSet_fastaData(sFastafile_open)
        Calculate_fasta_stat(dFastaDataSet)
        sFastafile_open.close()
        print("\nDone! it takes", round(time.time()-cTime, 4), "sec")

    elif len(args) == 2 and args[0] == 'SPLIT':
        sFastafile_open = open(args[1], 'rU')
        dFastaDataSet = Make_DataSet_fastaData(sFastafile_open)
        sys.exit(Split_single_to_multi(dFastaDataSet))

    else:
        print('''
==================================================================
Fasta & Fatsq modifier (ver. 0.2)

FQ2FA : Convert a FastQ file to a FastA file (stdout)
MF2SF : Split a multi-fasta to multiple single fasta files
STAT  : Calculate basic stats (#seq, bps, N50, %ATGCN) from a multi-fasta file
SPLIT : Split FastA to 60-len multi-lines FastA (stdout)

USAGE: fa_tool.py [FUNC] .fa or .fq

Aug, 2014 June
==================================================================''')
        quit()


if __name__ == "__main__":
    main(sys.argv[1:])

