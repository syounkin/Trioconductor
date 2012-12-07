setMethod("initialize", signature(.Object="TdtSet"),
	  function(.Object,
		   assayData=assayDataNew(calls=calls, ...),
		   calls=array(NA, dim=c(0,0,3)),
		   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
		   protocolData=phenoData[, integer(0)],
		   annotation=character(),
		   featureData=GenomeAnnotatedDataFrameFrom(assayData, annotation, genome=genome),
		   experimentData=new("MIAME"),
		   pedigree=Pedigree(),
		   genome=c("hg19", "hg18"), ...){
		  .Object <- callNextMethod(.Object,
					    assayData=assayData,
					    phenoData=phenoData,
					    featureData=featureData,
					    experimentData=experimentData,
					    annotation=annotation,
					    protocolData=protocolData,
					    pedigree=pedigree)
		  .Object@genome <- match.arg(genome)
		  return(.Object)
	  })

TdtSet <- function(pedigreeData=Pedigree(),
		   calls=array(NA, dim=c(0,0,3)), ...){
	if(ncol(calls) > 0) ids <- colnames(calls)
	if(!is.null(pedigreeData) & dim(calls)[[2]] > 0){
		if(!all(offspringNames(pedigreeData) %in% ids))
			pedigreeData <- pedigreeData[offspringNames(pedigreeData) %in% ids, ]
		fi <- match(fatherNames(pedigreeData), ids)
		mi <- match(motherNames(pedigreeData), ids)
		oi <- match(offspringNames(pedigreeData), ids)

		np <- nrow(trios(pedigreeData))
		nr <- nrow(calls)
		if(is(calls, "matrix")){
			gtArray <- initializeBigArray("baf", dim=c(nr, np, 3),
						      vmode="integer")
			dimnames(gtArray) <- list(rownames(calls),
						  offspringNames(pedigreeData),
						  c("F", "M", "O"))
			gtArray[, , 1] <- calls[, fi]
			gtArray[, , 2] <- calls[, mi]
			gtArray[, , 3] <- calls[, oi]
		} else {
			if(!is(calls, "array")) stop("calls must be a matrix or an array")
			gtArray <- calls
		}
	}
	new("TdtSet", pedigree=pedigreeData, calls=gtArray, ...)
}

setMethod("calls", "TdtSet", function(object) assayDataElement(object, "calls"))
