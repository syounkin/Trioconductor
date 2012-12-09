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
  index.vec <- which( colSums( apply( trios(ped), 1, FUN = function( row, id = id.vec ){ row %in% id } ) ) == 3 )
  ped.complete <- Pedigree( pedigreeInfo = trios(ped)[index.vec,] )
  return(ped.complete)
}

genoMat <- function( ped, data.mat ){
  ped.complete <- completeTrios( ped, as.character(rownames(data.mat)) )
  matrix.off <- data.mat[offspringNames(ped.complete),]
  matrix.fa <- data.mat[fatherNames(ped.complete),]
  matrix.ma <- data.mat[motherNames(ped.complete),]  
  matrix.trio <- weaveMat( mat.fa = matrix.fa, mat.ma = matrix.ma, mat.off = matrix.off ) 
  colnames(matrix.trio) <- c(t(cbind(fatherNames(ped.complete),motherNames(ped.complete),offspringNames(ped.complete))))
  rownames(matrix.trio) <- colnames(data.mat)
  return(matrix.trio)
}

weaveMat <- function( mat.fa, mat.ma, mat.off ){
  trio.mat <- matrix(c(rbind( mat.fa, mat.ma, mat.off )), nrow = ncol(mat.off), ncol = 3*nrow(mat.off), byrow = FALSE)
  colnames(trio.mat) <- c(t(cbind( rownames(mat.fa),rownames(mat.ma),rownames(mat.off))))
  rownames(trio.mat) <- colnames(mat.off)
 return( trio.mat )
}

geno.mat.fn <- function( ts, type = "holger" ){
  geno.array <- geno(ts)
  if( type == "holger" ){
    geno.out <- weaveMat( t(geno.array[,,"F"]),t(geno.array[,,"M"]),t(geno.array[,,"O"]) )
    return( t(geno.out) )
  }else{
    stop( "What format type would youlike for output?" )
  }
}
