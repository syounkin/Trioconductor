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

setMethod("show", "Pedigree", function(object){ print(head(object@trios) ) } )

setMethod("trios", signature(object="Pedigree"),
	  function(object) object@trios)

