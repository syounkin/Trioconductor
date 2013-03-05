setClass("PedClass", contains = "DataFrame",
         prototype = DataFrame(
           famid=paste0("fam",rep(1:3,3)),
           id=paste0("pseudo",1:9),
           fid=paste0("pseudo",1:9),
           mid=paste0("pseudo",1:9),
           sex=paste0("pseudo",1:9),
           dx=paste0("pseudo",1:9)
           ))

setClass("SNPTrioExperiment", contains="SummarizedExperiment", representation(pedigree="PedClass") )



## setOldClass("ff_array")
## setOldClass("ff_matrix")
## setClass("LogRratioSet", contains="eSet")
## #setClassUnion("matrixOrNULL", c("matrix", "NULL"))
## setClassUnion("arrayORff_array", c("array", "ff_array"))

## #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## #~~~~ TrioSet Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## setClass("tSet", contains="gSet",
## 	 representation(pedigree="Pedigree", "VIRTUAL"))

## setClass("TdtSet", contains="tSet")


## setClass("gTrio", contains="gSet",
##  	 representation(fatherPhenoData="AnnotatedDataFrame",
##  			motherPhenoData="AnnotatedDataFrame",
##  			pedigree="Pedigree", map="data.frame" ))

## setClass("iTrio", contains="gSet",
##  	 representation(fatherPhenoData="AnnotatedDataFrame",
##  			motherPhenoData="AnnotatedDataFrame",
##  			pedigree="Pedigree" ))

## #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## #~~~~ TrioSetList Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  setClass("TrioSetList", contains="gSetList",
##  	 representation(pedigree="Pedigree",
##  			fatherPhenoData="AnnotatedDataFrame",
##  			motherPhenoData="AnnotatedDataFrame"))
