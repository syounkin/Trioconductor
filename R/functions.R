aTDT <- function(geno, correct=FALSE){
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
