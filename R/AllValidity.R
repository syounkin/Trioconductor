validPedigree <- function(object){
	msg <- NULL
	if(nrow(trios(object)) > 0){
		if(!identical(colnames(trios(object)), c("F", "M", "O"))){
			msg <- "column names should be 'F', 'M', and 'O'"
			return(msg)
		}
		if(any(duplicated(offspringNames(object)))){
			msg <- "offspring identifiers must uniquely identify a trio"
			return(msg)
		}
		if(any(is.na(unlist(trios(object))))){
			msg <- "Missing values not allowed in pedigree"
			return(msg)
		}
		if(!all(c(is(fatherNames(object), "character"),
			  is(motherNames(object), "character"),
			  is(sampleNames(object), "character")))){
			msg <- "sample identifiers must be character strings (e.g., not factors)"
			return(msg)
		}
		if(!all(originalNames(allNames(object)) %in% originalNames(unlist(trios(object))))){
			msg <- "all 'individualId' in slot pedigreeIndex must correspond to an id in the trio slot"
			return(msg)
		}
		if(any(fatherNames(object) == motherNames(object))){
			msg <- "fatherNames can not be the same as the motherNames"
			return(msg)
		}
		if(any(fatherNames(object) == sampleNames(object))){
			msg <- "fatherNames can not be the same as the offspringNames"
			return(msg)
		}
		if(any(motherNames(object) == sampleNames(object))){
			msg <- "motherNames can not be the same as the offspringNames"
			return(msg)
		}
	}
	return(msg)
}

setValidity("Pedigree", function(object){
	msg <- validPedigree(object)
	if(is.null(msg)) return(TRUE) else return(msg)
})
