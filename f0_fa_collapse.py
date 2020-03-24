#!/usr/bin/env python3


import sys
import time


arg = sys.argv[1:]

############################################
###### Funtions

def CountSeqAbunt(sFile_ls):
    '''
    Count Seq Abunduncy from a Multifasta file.
    sFile_ls = list of .fa, striped, no blank line
    sLine_
    '''

    snSeqCount_dict = {}


    for sLine in sFile_ls:
        
        if not sLine.startswith(">"): #select seq data only
            sLine = sLine.upper()     #make all seq to upper()
            try:
                snSeqCount_dict[sLine] += 1 #If the seq is presence as a key
            except:                         # add value +1
                try:
                    sLineC_str = ''         #Reverse complement the seq
                    for sN in sLine:
                        if sN == "A":
                            sLineC_str += "T"
                        elif sN == "T":
                            sLineC_str += "A"
                        elif sN == "C":
                            sLineC_str += "G"
                        elif sN == "G":
                            sLineC_str += "C"
                        else:
                            sLineC_str += sN
                        #end of if
                    #end of for
                    sLineRC_str = sLineC_str[::-1]
                    snSeqCount_dict[sLineRC_str] += 1 #If the RC seq is presence as a key
                except:                              #add value +1
                    snSeqCount_dict[sLine] = 1       # either FS seq or RC seq is not presence as a key
                #end of for                          # add the seq as a new key
            #end of try
        #end of if
                    
    return snSeqCount_dict


def WriteOut(sFile_ls, snSeqCount_dict, OutFile):
    '''
    Write out the calcuated data from 'CountSeqAbunt'to a new fasta file.
    the out form as

    >Seq00001[two spaces][abundancy#]
    ATGT......
    >Seq00002[two spaces][abundancy#]
    GTCGA........
    '''

    nLine=len(sFile_ls)                     #Calculate the line number
    #nLine = sum(1 for line in sFile)
    #print(nLine)
    nCipher = len(str(nLine))               #Calculate the cipher of line number (10000 -> 5
    nDic = len(snSeqCount_dict)             #Calculate items.number of dict just for Verbose
    form = ">Seq%0."+str(nCipher)+"d %s"    #For text format [>Seq00001[two spaces][abundancy#]]
    j,k = 1,1


    for Seq, Num in snSeqCount_dict.items():
        OutFile.write(form % (j, Num))     #Print name of fasta
        OutFile.write("\n")                #Since .write did not add end as "\n"
        OutFile.write(Seq+"\n")            #Print seq data
        j += 1
        k += 1
        if k/nDic >= 0.1:                  #Verbose for processing
            k = 0
            print(str(int(j/nDic*100))+"%")
        #end of if
    #end of for
#end of def

####################################################
####### Main


def main(arg):

    if len(arg) == 2:

        print("File to be collapesd -- %s" % arg[0])
        print("Collaped file out as -- %s" % arg[1])

        print("start load %s" % arg[0])
        sFile = open(arg[0], "rU")
        sFile_ls = sFile.readlines()                          #file to list
        sFile_ls = [i.strip() for i in sFile_ls if i != "\n"] #eliminate blank line & strip()
        snSeqCount_dict = CountSeqAbunt(sFile_ls)
        print("end load %s, takes %s sec" % (arg[0], int(time.time()-ctime)))
        print("start output to %s" % arg[1])
        OutFile = open(arg[1], "w")
        WriteOut(sFile_ls, snSeqCount_dict, OutFile)
        print("end ountput to %s, takes %s sec" % (arg[1], int(time.time()-ctime)))
        OutFile.close()

    else:
        print('''
============================================================================
Fasta Collapser (ver 0.1)

Usage : python3 fa_collapse.py [input.fa] [output.fa]

Aug, 2014 June
============================================================================''')
        quit()
    #end of if
#end of main

     
if __name__ =="__main__":
    ctime = time.time()
    main(arg)

