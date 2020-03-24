#!/usr/bin/env python3


import glob
import sys

################################################################
########## args

arg = sys.argv[1:]

if len(arg) == 0:
    print()
    print('''
    Multiple VCF file manager v0.1
    This tool could merge, compare & parse the Variant Call Format (VCF) files.

    merge   : Merging multiple VCF files to single output, contain only common SNP info in all samples
              Input should be provided as regular expression (e.g. [7-9]*.vcf)
              output will be [STDOUT]
    compare : Compare two VCF files (also 'merge' output files) & make 3 output files containing
              common, sample1-specific, & sample2-specific SNPs for each
    parse   : Parse SNP data from multiple vcf files to a reference vcf file (also 'merge' output files).
              Inputs should be provided as regular expression
              output will be [STDOUT]

    Usage
    python3 vcf.parse.py merge [inputs]
    python3 vcf.parse.py compare [input1] [input2]
    python3 vcf.parse.py parse [ref. vcf] [inputs]
           ''')
    sys.exit()
else:
    pass
#end of if

################################################################
######### Methods

def VcfMerge(in1_ls):
    """
    Multiple .vcf merge to single .vcf files
    Call multiple .vcf files with glob;
    Call snp site and SNP info from multiple files,
    and only commonly exist and the sam SNPs in all .vcf files are printed
    "len(n) = num" ; len(set(n(1:))) == 1

    Oct, 2014 by June
    """

    num = len(in1_ls)
    k = 1
    dic = {}

    for i in in1_ls:
            fi = open(i)
            for j in fi:
                    j = j.strip().split('\t') if not j.startswith('##') else ''
                    if j != '' and k == 1:
                            site = j[0]+'@'+str(j[1])
                            dic[site] = [j[3], j[4]]
                    elif j != '':
                            site = j[0]+'@'+str(j[1])
                            snp = j[4]
                            try:
                                    dic[site].append(snp)
                            except KeyError:
                                    pass
                            #end of try
                    #end of if
            k += 1
            #end of for
    #end of for

    print("#contg", "site", "ref", "snp", sep="\t")
    for l,n in dic.items():
            if len(n) == num and len(set(n[1:])) == 1:
                    site = l.split('@')
                    print(site[0], site[1], sep="\t", end="\t")
                    for m in n:
                            print(m, end="\t")
                    #end of for
                    print()
            #end of if
    #end of for

def FindCommSNP(vcf1, vcf2):
    """
    Comparing tools of two merged .vcf files (vcf1, vcf2)
    Make 3 outputs;
    Common SNPs (out1), specific SNPs for each .vcf file (out2 & out3)
    will be calculated

    Oct. 2014 June
    """
    dic1 = {}
    dic2 = {}
    out1 = open('common.snp.vcf', 'w')
    out2 = open('vcf1.only.snp.vcf', 'w')
    out3 = open('vcf2.only.snp.vcf', 'w')

    for i in vcf1:
        i = i.strip().split('\t')
        site = i[0]+"@"+i[1]
        dic1[site] = i[2:]
    #end of for
    for j in vcf2:
        j = j.strip().split('\t')
        site = j[0]+"@"+j[1]
        dic2[site] = j[2:]
    #end of for

    print("#contg", "site", "ref", "snp", sep="\t", file=out1)
    print("#contg", "site", "ref", "snp", sep="\t", file=out2)

    common = []
    for l,m in dic1.items():
        try:
            if dic1[l][-1] == dic2[l][-1]:
                key = l.split('@')
                print(key[0],key[1],dic1[l][0], dic1[l][1], sep="\t", file=out1)
                common.append(l)
            else: pass
            #end of if
        except KeyError:
            pass
        #end of try
    #end of for
    for z in common:
        del dic1[z]
        del dic2[z]
    #end of for
    for l, m in dic1.items():
        l = l.split('@')
        print(l[0], l[1], file=out2, sep="\t", end="\t")
        for x in m:
            print(x, file=out2, end="\t")
        print(file=out2)
        #end of for
    #end of for
    out1.close()
    out2.close()

    print("#contg", "site", "ref", "snp", sep="\t", file=out3)
    for n in dic2.keys():
        key = n.split('@')
        print(key[0], key[1], file=out3, sep="\t", end="\t")
        for y in dic2[n]:
            print(y, file=out3, end="\t")
        print(file=out3)
        #end of for
    #end of for
    out3.close()

def SnpTargetParse(vcf3, in2_ls):
    """
    .vcf parser to find SNP variations from multiple .vcf files
    work like 'join' function in Shell

    Call location, reference base, SNP base from 'fi1' opened file,
    which offer location info to call SNP from other files

    Call location, SNP info from fi2 file lists, & add the SNP info to associataed location of 'fi1'

    Output appears like below
    chr site ref SNP_fi1 SNP_fi2-1 SNP_fi2-2 ....

    Oct. 2014 by june
    """
    fi1 = open('6.work.vcf')
    fi2_ls = glob.glob('./3_vcf/[7-9]*.vcf')

    dic = {}
    for i in fi1:
        i = i.strip().split('\t')
        site = i[0]+":"+str(i[2])
        dic[site] = [i[3], i[4]]
    #end of for

    count = 3
    for fi2 in fi2_ls:
        fi2 = open(fi2)
        for l in fi2:
            if not l.startswith("##"):
                l = l.strip().split('\t')
                site2 = l[0][:-2]+":"+str(l[1])
                try:
                    dic[site2].append(l[4])
                except KeyError:
                    pass
                #end of try
            #end of if
        for m in dic.values():
            if len(m) == count-1:
                m.append(m[0])
            elif len(m) == count:
                pass
            else:
                print("Error")
                break
            #end of if
        #end of for
        count += 1
    #end of for

    print("#gene", "site", "ref", "ref_SNP", "SNP1-1", "...", sep="\t")
    for x, y in dic.items():
        x = x.split(':')
        print(x[0],x[1],sep="\t", end="\t")
        for i in y:
            print(i, end='\t')
        print()
        #end of for
    #end of for

#####################################################################
############### Main

def main(arg):
    if len(arg) == 2 and arg[0] == 'merge':
        in1_ls = glob.glob(arg[1])
        VcfMerge(in1_ls)

    elif len(arg) == 3 and arg[0] == 'compare':
        vcf1 = open(arg[1])
        vcf2 = open(arg[2])
        FindCommSNP(vcf1, vcf2)

    elif len(arg) == 3 and arg[0] == 'parse':
        vcf3 = open(arg[1])
        in2_ls = glob.glob(arg[2])
        SnpTargetParse(vcf3, in2_ls)

    else:
        print()
        print('''
        Multiple VCF file manager v0.1
        This tool could merge, compare & parse the Variant Call Format (VCF) files.

        merge   : Merging multiple VCF files to single output, contain only common SNP info in all samples
                  Input should be provided as regular expression (e.g. [7-9]*.vcf)
                  output will be [STDOUT]
        compare : Compare two VCF files (also 'merge' output files) & make 3 output files containing
                  common, sample1-specific, & sample2-specific SNPs for each
        parse   : Parse SNP data from multiple vcf files to a reference vcf file (also 'merge' output files).
                  Input should be provided as regular expression (e.g. [7-9]*.vcf)
                  output will be [STDOUT]

        Usage
        python3 vcf.parse.py merge [inputs]
        python3 vcf.parse.py compare [input1] [input2]
        python3 vcf.parse.py parse [ref.vcf] [inputs]
               ''')
    #end of if

if __name__ == '__main__':
    main(arg)
