annotatedDataFrameFromArray <- function(object, byrow=FALSE, ...){
	if(dim(object)[[3]] > 0){
		object <- object[, , 1, drop=TRUE]
		res <- Biobase:::annotatedDataFrameFromMatrix(object, byrow=byrow, ...)
	} else res <- Biobase:::annotatedDataFrameFromMatrix(matrix(), byrow=byrow, ...)
	return(res)
}

GenomeAnnotatedDataFrameFromArray <- function(object, annotationPkg, genome, ...){
	## coerce to matrix
	dims <- dim(object)
	is.array <- length(dims) == 3
	if(is.array){
		res <- oligoClasses:::GenomeAnnotatedDataFrameFromMatrix(object[, , 1], annotationPkg=annotationPkg, genome=genome, ...)
	} else {
		##dim(object) <- dim(object)[c(1,2)]
		res <- oligoClasses:::GenomeAnnotatedDataFrameFromMatrix(object,
									 annotationPkg=annotationPkg,
									 genome=genome, ...)
	}
	res
}


setMethod("annotatedDataFrameFrom", signature(object="ff_array"),
	  annotatedDataFrameFromArray)

setMethod("annotatedDataFrameFrom", signature(object="array"),
	  annotatedDataFrameFromArray)

## this should probably be moved to VanillaICE, then imported by MinimumDistance
setMethod("GenomeAnnotatedDataFrameFrom", signature(object="array"),
	  function(object, annotationPkg, genome, ...){
		  GenomeAnnotatedDataFrameFromArray(object, annotationPkg=annotationPkg, genome=genome, ...)
	  })



