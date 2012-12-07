test_TdtSet_construction <- function(){
	## an empty TdtSet
	checkTrue(validObject(TdtSet()))

	path <- system.file("extdata", package="MinimumDistance")
	load(file.path(path, "logRratio.rda"))
	load(file.path(path, "pedigreeInfo.rda"))
	calls <- matrix(sample(1:3, 550*6, replace=TRUE), 550, 6)
	colnames(calls) <- as.character(unlist(pedigreeInfo))
	ped <- Pedigree(pedigreeInfo)
	tdt.obj <- TdtSet(pedigreeData=ped, calls=calls)

	## this is an array
	genotypes <- calls(tdt.obj)
	## father, mother, offspring genotypes for first 5 snps
	genotypes[1:5, , ]

	## need set validity method for this class

}
