aTDT.fn <- function(geno, correct=FALSE){
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
  return( list(stat=stat, pval=pval, transMinor=transMinor, transMajor=transMajor ))
}

trans.tab <- function( object ) {
  T <- sum(object[c("101","011","122","212","111","112","112")],na.rm=TRUE)
  U <- sum(object[c("100","010","121","211","111","110","110")],na.rm=TRUE)
  if( T+U >= 50 ){
    return(binom.test(x = T,n = T+U,p = 0.5, alternative = "greater")$p.value)
  }else{
    return(NA)
  }
}

trans.rate <- function( object ) {
  T <- sum(object[c("101","011","122","212","111","112","112")],na.rm=TRUE)
  U <- sum(object[c("100","010","121","211","111","110","110")],na.rm=TRUE)
  if( T+U >= 50 ){
    return(binom.test(x = T,n = T+U,p = 0.5))
  }else{
    return(NA)
  }
}

get.est <- function(obj){
  if(!is.na(obj)){
    return(obj$estimate)
  }else{
    return(NA)
  }
}

get.ci <- function(obj){
  if(!is.na(obj)){
    return(obj$conf.int)
  }else{
    return(NA)
  }
}

make.files.for.cpp <- function(obj, fileroot){
  ped.df <- as(pedigree(obj), "data.frame")
  geno.mat <- as(as(geno(obj),"numeric"), "matrix")
  rownames(ped.df) <- ped.df$id
  pedfile <- cbind(ped.df[rownames(geno.mat),], geno.mat)[c(t(as.matrix(completeTrios(obj)))),]
  write.table(pedfile, file = paste0(fileroot,".ped"), quote = FALSE, col.names = FALSE, row.names = FALSE)
    write.table(start(rowData(obj)), file = paste0(fileroot,".map"), quote = FALSE, col.names = FALSE, row.names = FALSE)
  nsnps <- length(rowData(obj))
  write.table(rep(1,nsnps), file = paste0(fileroot,".weights"), quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(data.frame(name="pseudoblock", start=1, end=nsnps), file = paste0(fileroot,".blocks"), quote = FALSE, col.names = FALSE, row.names = FALSE)
  }

CountTU <- function( x ){
  T <- sum(c(x["011"],x["101"],x["111"],x["112"],x["112"],x["122"],x["212"]), na.rm = TRUE )
  U <- sum(c(x["010"],x["100"],x["111"],x["110"],x["110"],x["121"],x["211"]), na.rm = TRUE )
  return(c(T,U))
}

TU.fish <- function( TU.vec ){
  if(!identical(c("T.case", "U.case", "T.con", "U.con"),names(TU.vec))) return( "TU vector must have the correct names.")
  fish <- fisher.test(matrix(TU.vec, nrow = 2, ncol = 2, byrow = TRUE), alternative = "greater")
  return(fish)
}


