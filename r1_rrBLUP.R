## rrBLUP progress log
## The input file ends with "BLUP.geno" is a output of vF_vcf2blup.py from vcf file
## Calculate all phenotypes provided and outputs "ans1.blup.txt"
## the output file can be manipulted by vB_plink2qqman.py to QQMAN-convinient files.
## Mar 27th 2016
## Junesk9


library(rrBLUP)
geno <- read.table("AT.E180plus.2nd.SNP.BLUP.geno",header=T,as.is=T,check.name=FALSE)

m <- ncol(geno)
geno2 <- geno[,4:m]
A <- A.mat(t(geno2),n.core=30,impute.method="EM")

pheno <- read.table("99_flower.pheno.txt",header=T,as.is=T,check.name=FALSE)
pheno2 <- pheno[,-2]
colnames(pheno2)[1] = "ecoid"
eco <- colnames(geno)[c(-1:-3)]
pheno3 <- subset(pheno2, ecoid %in% eco)
pheno3[pheno3=="-9"] <- NA

ans2 <- GWAS(pheno=pheno3, geno=geno, P3D=TRUE, n.core=35,K=A, plot=FALSE)
m <- ncol(ans2)
for (i in 4:m){
	ans2$newcol <- p.adjust(10^-ans2[,i],method="fdr") # temporary, add FDR info to a new colomn "newcol"
	colnames(ans2)[ncol(ans2)] <- paste("fdr",i-3,sep="") # change the name of new column "newcol" > fdrXX
}

write.table(ans2,"ans1.blup.out.txt",sep="\t",row.names=F,quote=F) 



