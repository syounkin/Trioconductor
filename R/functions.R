Pedigree <- function(pedigreeInfo,
		     fatherIds=character(),
		     motherIds=character(),
		     offspringIds=character()){
	if(!missing(pedigreeInfo)){
          if( !all( c("F","M","O") %in% names(pedigreeInfo) ) ){
            stop("pedigreeInfo must contain columns named F, M, and O for father, mother and offspring ids.") }else{
		msg <- "pedigreeInfo must be a data.frame with column names 'F', 'M', and 'O'"
		if(!is(pedigreeInfo, "data.frame"))
			stop(msg)
		trios <- data.frame(F=as.character(pedigreeInfo$F),
				    M=as.character(pedigreeInfo$M),
				    O=as.character(pedigreeInfo$O),
				    stringsAsFactors=FALSE)
		allIds <- as.character(unlist(trios))
              }
	} else {
		fatherIds <- as.character(fatherIds)
		motherIds <- as.character(motherIds)
		offspringIds <- as.character(offspringIds)
		trios <- data.frame(F=fatherIds,
				    M=motherIds,
				    O=offspringIds,
				    stringsAsFactors=FALSE)
		allIds <- c(fatherIds, motherIds, offspringIds)
	}
	trio.index <- as.integer(matrix(seq_len(nrow(trios)), nrow(trios), 3, byrow=FALSE))
	memberId <- rep(c("F", "M", "O"), each=nrow(trios))
	pedigreeIndex <- data.frame(individualId=allIds,
				    memberId=memberId,
				    index.in.pedigree=trio.index,
			    stringsAsFactors=FALSE)
	rownames(pedigreeIndex) <- NULL
	new("Pedigree", trios=trios, trioIndex=pedigreeIndex)
}
