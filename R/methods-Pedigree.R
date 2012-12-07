setMethod("initialize", signature(.Object="Pedigree"),
	  function(.Object,
		   trios=data.frame(F=character(),
		   M=character(),
		   O=character(), stringsAsFactors=FALSE),
		   trioIndex=data.frame(individualId=character(),
		   memberId=character(), index.in.pedigree=integer(), stringsAsFactors=FALSE),
		   ...){
		  callNextMethod(.Object, trios=trios, trioIndex=trioIndex, ...)
	  })

setMethod("show", signature(object="Pedigree"), function(object){
  cat(paste0("This pedigree object contains ", nrow(trios(object)), " complete trios.\nFor access to the data frame use the trios() accessor function.\n"))
} )
#~~~~~~~~~~~~~~~~
setMethod("trios", signature(object="Pedigree"), function(object) object@trios)
setMethod("trioIndex", signature(object="Pedigree"), function(object) object@trioIndex)
#~~~~~~~~~~~~~~~~
setMethod("offspringNames", signature(object="Pedigree"), function(object) trios(object)$O)
setMethod("fatherNames", signature(object="Pedigree"), function(object) trios(object)$F)
setMethod("motherNames", signature(object="Pedigree"), function(object) trios(object)$M)
setMethod("sampleNames", signature(object="Pedigree"), function(object) offspringNames(object))
setMethod("allNames", signature(object="Pedigree"), function(object) unique(trioIndex(object)$individualId))
setMethod("completeTrios", signature(object = "Pedigree", id.vec = "character"), function( object, id.vec ){ completeTrios.fn(object,id.vec) } )
