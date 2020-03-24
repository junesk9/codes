Toshimori Tools
(Last modified: Nov.18th.2014)

############ FA or FQ file handling

f0. fa_collapse.py
	From NGS-based short-read .fa file,
	collapse the same sequence in .fa file & count the same seqs
	(also consider reverse-complement)

	For, sRNA analysis; though not very useful, since many sRNA analysis tools embeding this function.

f1. fa_modi.py
	Harboring many basic tools for .fa .fq file handling
	FQ2FA : Convert .fq file to .fa
	MF2SF : Convert a multi-fasta file to multiple single fasta files
	STAT : Calculate basic stat of .fa
	SPLIT : fold single line .fa to specific length (default: 60); [fold] in bash is more useful

f2. fa_select.py
	Collect specified sequences from a multiple .fa file
	
f3. fq_count.py
	Caculate total bp. and total seq number from a .fq file

f4. chr_info_extr.py
	From a multi-fa file,
	extract simplified info as [name len]
	
	For, bedtools, vcftools, et cetra.

############ sRNA data manipulation

s1. fa_size_collect.py
	From NGS-based short read .fa file,
  	only seq has given size retained (default: 24)

s2. make.mirGFF.py
	parse ShortStack (mirRNA prediction suite) output to miRNA info as .gff format
	each miRNA have three lines. one for premature (stem-loop structure), another for mature-5p, the other for mature-3p

s3. eliminate_miR_fa.py
	miR sequence eliminator
	Eliminate mature miR sequence from a multi-fasta file
	output as STDOUT

s4. SS_output_parser.py
	Parse ShortStack output 'MIRNAs' folder
	to single output txt file
	
	excute at the /ShortStack_***********/

############# VCF data manipulation

v1. vcf_tools.py
	Useful tools for .vcf (Variant Call File) data manipulation
	
	merge   : Merging multiple VCF files to single output, contain only common SNP info in all samples
                  Input should be provided as regular expression (e.g. [7-9]*.vcf)
                  output will be [STDOUT]
        compare : Compare two VCF files (also 'merge' output files) & make 3 output files containing
                  common, sample1-specific, & sample2-specific SNPs for each
        parse   : Parse SNP data from multiple vcf files to a reference vcf file (also 'merge' output files).
                  Inputs should be provided as regular expression
                  output will be [STDOUT]

############# BSSEQ data manipulation

m1. wig2input.py
	Some published  BSSEQ data distributed by wiggle format.
	since this tools convert .wig to .input format for 5mC.comparing pipeline

m2. bwa-bed.bias-survey.py
	For bwa-meth output bedGraph file.
	surveying frequency of read depth on each Cytosine output
	
	Useful to check the average read depth of mapping result
	& decise cut-off depth for reliable cytosines

m3. 5mC_genome_profile.py
	Survey average 5mC rate on chromosome with specific window size
	Use .input file prepared for 5mC comparing pipeline

m4. find.comm.DMR.py
	Comparing two DMR files as ouput of 5mC comparing pipelines &
	output common DMRs with each coverage & same DMR context

m5. find.DMR.prox.gene.py
	Input as DMR file from 5mC comparing pipeline or [find.comm.DMR] & .gff files
	call genes in proximal range (default: 2k) around a DMR
	.gff files would be manipulated by purpose

m6. call.DMR.5mc.mc.py  **Take long time (more than hours)
	Input as DMR file from 5mC comparing pipelines or [find.comm.DMR] & multiple .input files for 5mC comparing pipelines
	call 5mC rate at refered DMR regions, single output have three column by seq context (CG, CHG & CHH)
	
	Equiped multiprocessing; the #num of threads as #chr is needed.

m7. parse.input.4boxplot.py
	multiple inputs from [call.DMR.5mc.mc.py] are parsed into three outpus by seq context,
	make easy to adapt R.boxplot

m8. dmr.chr.dist.py
	Survey Chromosome distributions of # of assigned both DMC & DMR by sequence context.
	make easy to adapt R.matplot
	Not very ubiquatous; Now only adaptable for A. thaliana
	input as DMC&DMR from DMC&DMR analysis, common DMC from m9, common DMR from m4

m9. find.comm.dmc.py
	Find common DMC from TWO DMC files with TWO criteria
	1. same chr, sequence
	2. the 5mC rates are both positive or negative numbers

mA. input2wig.py
	Convert a 5mC input file to three WIGGLE file by seq context (CG, CHG, CHH)
	Nov.2014 now, only A. thaliana data is adaptable

mB. bed2input.py
	Convert a bwa-meth output bed file to a 5mC input file
	Attend coutoff of read depth (dafault = 4)
