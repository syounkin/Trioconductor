### R code from vignette source './trioClasses.Rnw'

###################################################
### code chunk number 1: trioClasses.Rnw:3-4
###################################################
  options(width=70, continue = " ")


###################################################
### code chunk number 2: packages
###################################################
library("trioClasses")
library("vcf2R")
data(targets)
data("BMP4-european-all.sm")
df <- DataFrame( dx = rep(1,10), sex = rep(1,10), row.names = paste0("sub",1:10) )
class(df)
class(sm)
logR <- matrix( runif(130), nrow = 13, ncol = 10 )
baf <-  matrix( runif(130), nrow = 13, ncol = 10 )


###################################################
### code chunk number 3: peddf
###################################################
ped.df <- data.frame( pedid = c(1,1,1,2,2,2,3,3,3,4), id = paste0("sub",1:10), fid = c(0,0,"sub1",0,0,"sub4",0,0,"sub7",0), mid = c(0,0,"sub2",0,0,"sub5",0,0,"sub8",0), sex = c(1,2,2,1,2,2,1,2,2,1), dx = 1, stringsAsFactors=FALSE )


###################################################
### code chunk number 4: snpEx
###################################################
snpEx <- new("SNPTrioExperiment", pedigree = new("PedClass",ped.df), assays = SimpleList(geno = sm[1:13,1:10], logR = logR, baf = baf), colData = df, rowData = targets.gr ) 


###################################################
### code chunk number 5: classes
###################################################
class(snpEx@pedigree)
class(geno(snpEx))
class(logR(snpEx))
class(baf(snpEx))


###################################################
### code chunk number 6: methods
###################################################
pedigree(snpEx)
completeTrios(snpEx)


