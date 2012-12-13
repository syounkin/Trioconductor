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
  matrix.trio <- weaveMat( matrix.fa, matrix.ma, matrix.off ) 
  colnames(matrix.trio) <- c(rbind(fatherNames(ped.complete),motherNames(ped.complete),offspringNames(ped.complete)))
  rownames(matrix.trio) <- colnames(data.mat)
  return(matrix.trio)
}

weaveMat <- function( x, y, z ){
trio.mat <- matrix(c(t(cbind(x,y,z))),nrow = ncol(x), ncol = 3*nrow(x))
#trio.mat <- matrix(c(cbind( mat.fa, mat.ma, mat.off )), nrow = ncol(mat.off), ncol = 3*nrow(mat.off), byrow = TRUE)
  #colnames(trio.mat) <- c(t(cbind( rownames(mat.fa),rownames(mat.ma),rownames(mat.off))))
  #rownames(trio.mat) <- colnames(mat.off)
 return( trio.mat )
}

geno.mat.fn <- function( ts, type = "holger" ){
  geno.array <- geno(ts)
  id <- trios(gTrio.obj)$O
  fid <- trios(gTrio.obj)$F
  mid <- trios(gTrio.obj)$M
  names(fid) <- id
  names(mid) <- id
 
  if( type == "holger" ){
    geno.out <- weaveMat( t(geno.array[,,"F"]),t(geno.array[,,"M"]),t(geno.array[,,"O"]) )
    colnames(geno.out) <- c(t(cbind(fid,mid,id)))
    rownames(geno.out) <- rownames(geno.array[,,"F"])
    return( t(geno.out) )
  }else{
    stop( "What format type would you like for output?" )
  }
}

make.unique2 <- function(names, sep="___DUP") make.unique(names, sep)

originalNames <- function(names){
	if(length(names) ==0) return(names)
	sep <- formals(make.unique2)[["sep"]]
	index <- grep(sep, names)
	if(length(index) > 0) names[index] <- sapply(names[index], function(x) strsplit(x, sep)[[1]][[1]])
	return( names )
}

# Written by Holger Schwender, originally named aTDTchunk
aTDT <- function(gTrio, correct = FALSE){
        geno <- getGeno( gTrio, type = "holger" )
	n.row <- nrow(geno)
	dad <- geno[seq.int(1, n.row, 3),, drop=FALSE]
	mom <- geno[seq.int(2, n.row, 3),, drop=FALSE]
	kid <- geno[seq.int(3, n.row, 3),, drop=FALSE]
	dad <- 100 * dad + 10 * mom + kid
	tmp1 <- colSums(dad==11 | dad==101, na.rm=TRUE)
	tmp2 <- colSums(dad==122 | dad==212, na.rm=TRUE)
	tmp3 <- colSums(dad==111, na.rm=TRUE)
	tmp4 <- 2 * colSums(dad==112, na.rm=TRUE)
	transMinor <- tmp1 + tmp2 + tmp3 + tmp4
	tmp1 <- colSums(dad==10 | dad==100, na.rm=TRUE)
	tmp2 <- colSums(dad==121 | dad==211, na.rm=TRUE)
	tmp4 <- 2 * colSums(dad==110, na.rm=TRUE)
	transMajor <- tmp1 + tmp2 + tmp3 + tmp4
	tmp1 <- transMinor - transMajor
	if(correct)
		tmp1 <- abs(tmp1) - 1
	stat <- tmp1 * tmp1 / (transMinor + transMajor)
      	pval <- pchisq(stat, 1, lower.tail=FALSE)
        return( data.frame( stat=stat, pval=pval, transMinor=transMinor, transMajor=transMajor ) )
}
