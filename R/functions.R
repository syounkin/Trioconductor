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

genoMat <- function( ped, data.mat ){
  ped.complete <- completeTrios( ped, rownames(data.mat) )
  matrix.off <- data.mat[offspringNames(ped.complete),]
  matrix.fa <- data.mat[fatherNames(ped.complete),]
  matrix.ma <- data.mat[motherNames(ped.complete),]  
  matrix.trio <- t(matrix(c(rbind( matrix.fa, matrix.ma, matrix.off )), nrow = 3*nrow(matrix.off), byrow = TRUE))
  colnames(matrix.trio) <- c(t(cbind(fatherNames(ped.complete),motherNames(ped.complete),offspringNames(ped.complete))))
  rownames(matrix.trio) <- colnames(data.mat)
  return(matrix.trio)
}
