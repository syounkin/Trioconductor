make.unique2 <- function(names, sep="___DUP") make.unique(names, sep)
originalNames <- function(names){
	if(length(names) ==0) return(names)
	sep <- formals(make.unique2)[["sep"]]
	index <- grep(sep, names)
	if(length(index) > 0) names[index] <- sapply(names[index], function(x) strsplit(x, sep)[[1]][[1]])
	names
}

completeTrios.fn <- function( ped, id.vec ){
  if( !is.character(id.vec) ) stop( "ID vector must be character class.")
  index.vec <- which(colSums(apply( trios(ped), 1, FUN = function( row, id = id.vec ){row %in% id }))==3)
  ped.complete <- Pedigree( pedigreeInfo = trios(ped)[index.vec,] )
  return(ped.complete)
}

