setClass("Pedigree", representation(trios="data.frame", trioIndex="data.frame"))

setOldClass("ff_array")
setOldClass("ff_matrix")
setClass("LogRratioSet", contains="eSet")
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("arrayORff_array", c("array", "ff_array"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~ TrioSet Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setClass("tSet", contains="gSet",
	 representation(pedigree="Pedigree", "VIRTUAL"))

setClass("TdtSet", contains="tSet")


setClass("TrioSet", contains="gSet",
 	 representation(fatherPhenoData="AnnotatedDataFrame",
 			motherPhenoData="AnnotatedDataFrame",
 			pedigree="Pedigree",
 			mindist="matrixOrNULL"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~ TrioSetList Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 setClass("TrioSetList", contains="gSetList",
 	 representation(pedigree="Pedigree",
 			fatherPhenoData="AnnotatedDataFrame",
 			motherPhenoData="AnnotatedDataFrame"))
