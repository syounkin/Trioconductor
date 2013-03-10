setMethod("initialize", "SNPTrioExperiment", function(.Object, pedigree, ... ){
  .Object@pedigree <- pedigree
  .Object <- callNextMethod()
  .Object
})

setMethod("SNPTrioExperiment", signature("SummarizedExperiment", "PedClass"), function(se, pedigree){
  new("SNPTrioExperiment", pedigree, se)
})

setMethod("show", signature(object="SNPTrioExperiment"), function(object){
  callNextMethod()
  cat("pedigree(", nrow(pedigree(object)), "): famid id fid mid sex dx\n", sep = "")
})

setMethod("geno", signature(object="SNPTrioExperiment"), function(object) {
  geno <- t(assays(object)$geno)
  colnames(geno) <- as.character(names(rowData(object)))
  rownames(geno) <- as.character(rownames(colData(object)))
  geno
})

setMethod("logR", signature(object="SNPTrioExperiment"), function(object) assays(object)$logR )

setMethod("baf",  signature(object="SNPTrioExperiment"), function(object) assays(object)$baf )

setMethod("pedigree",  signature(object="SNPTrioExperiment"), function(object) object@pedigree )

setMethod("completeTrios",  signature(object="SNPTrioExperiment"), function(object){
  ped <- pedigree(object)
  trios <- trios(ped)
  ids <- colnames(object)
  index <- (trios[,"id"] %in% ids) & (trios[,"fid"] %in% ids )& (trios[,"mid"] %in% ids)
  trios[index,]
})

setMethod("GenoTrio",  signature(object="SNPTrioExperiment"), function(object){
  ct <- completeTrios(object)
  geno <- geno(object)
  list( O = geno[as.character(ct[,"id"]),], F = geno[as.character(ct[,"fid"]),], M = geno[as.character(ct[,"mid"]),] )
})

setMethod("ctcbind", signature( object = "list"), function( object ){
  holger <- with( object, matrix(c(t(cbind( as(F,"matrix"), as(M,"matrix"), as(O,"matrix")))), byrow = TRUE, nrow = 3*nrow(O), ncol = ncol(O)))
  x <- ifelse( holger == 01, 0L, ifelse( holger == 02, 1L, ifelse( holger == 03, 2L, NA )))
  x
})


setMethod("[", "SNPTrioExperiment", function( x, i, j, ..., drop = TRUE ){
  se <- as( x, "SummarizedExperiment" )
  if( !missing(i) & !missing(j)){
    se <- se[i,j]
  }else if( !missing(i) & missing(j)){
    se <- se[i,]
  }else if( missing(i) & !missing(j)){
    se <- se[,j]
  }else {
    se <- se [,]
  }
  return(SNPTrioExperiment(se, pedigree(x)))
})

setMethod("parents",  signature(object="SNPTrioExperiment"), function(object){
  with(as(pedigree(object),"data.frame"), unique(c(as.character(fid), as.character(mid))))
})

setMethod("MAF", signature(object="SNPTrioExperiment"), function(object){
  with(col.summary(geno(ste[,which(colnames(ste) %in% parents(ste))])),MAF)
})

setMethod("aTDT", signature(object="SNPTrioExperiment"), function(object){
  geno <- as( object, "matrix" )
  aTDT.fn(geno)
})

setMethod("aTDT", signature(object="matrix"), function(object){
  geno <- object
  aTDT.fn(geno)
})

setAs( from = "SNPTrioExperiment", to = "matrix", function(from){
  gtrio <- GenoTrio(from)
  holger <- ctcbind(gtrio)
  return(holger)
})

## setMethod("TransCount", signature( object = "SNPTrioExperiment", region = "GRanges"), function( object, region ){
##   index <- subjectHits(findOverlaps(region, rowData(object)))
##   t <- numeric(6)
##   if( length(index) != 0 ){
##       ste <- object[index,]
##       gtrio <- GenoTrio(ste)
##       t[1] <-   with( gtrio, sum( F == 1 & M == 2 & O == 2, na.rm = TRUE) )
##       t[2] <-   with( gtrio, sum( F == 2 & M == 1 & O == 2, na.rm = TRUE) )
##       t[3] <-   with( gtrio, sum( F == 2 & M == 3 & O == 3, na.rm = TRUE) )
##       t[4] <-   with( gtrio, sum( F == 3 & M == 2 & O == 3, na.rm = TRUE) )
##       t[5] <-   with( gtrio, sum( F == 2 & M == 2 & O == 2, na.rm = TRUE) )
##       t[6] <- 2*with( gtrio, sum( F == 2 & M == 2 & O == 3, na.rm = TRUE) )
##       return( sum(t, na.rm = TRUE ))
##   }else{
##     return(NA)
##   }
## })

setMethod("TransCount", signature( object = "SNPTrioExperiment", region = "GRanges"), function( object, region ){
  minor <- matrix(0,nrow=length(region), ncol = 6)
  major <- matrix(0,nrow=length(region), ncol = 6)
  mendel <- matrix(0,nrow=length(region), ncol = 4)
  for( i in 1:length(region) ){
    gr.row <- region[i]
    ste <- object[subjectHits(findOverlaps(gr.row, rowData(object))),]
    gtrio <- GenoTrio(ste)
#    with( gtrio, {
      F <- as(with(gtrio,F), "numeric")
      M <- as(with(gtrio,M), "numeric")
      O <- as(with(gtrio,O), "numeric")
      
      mendel[i,1] <- sum( F == 0 & M == 0 & ( O == 1 | O == 2 ), na.rm = TRUE)
      mendel[i,2] <- sum( F == 2 & M == 2 & ( O == 0 | O == 1 ), na.rm = TRUE)
      mendel[i,3] <- sum( F == 0 & M == 2 & ( O == 0 | O == 2 ), na.rm = TRUE)
      mendel[i,4] <- sum( F == 2 & M == 0 & ( O == 0 | O == 2 ), na.rm = TRUE)      
    
      minor[i,1] <- sum( F == 0 & M == 1 & O == 1, na.rm = TRUE)
      major[i,1] <- sum( F == 0 & M == 1 & O == 0, na.rm = TRUE)

      minor[i,2] <- sum( F == 1 & M == 0 & O == 1, na.rm = TRUE)
      major[i,2] <- sum( F == 1 & M == 0 & O == 0, na.rm = TRUE)

      minor[i,3] <- sum( F == 1 & M == 2 & O == 2, na.rm = TRUE)
      major[i,3] <- sum( F == 1 & M == 2 & O == 1, na.rm = TRUE)

      minor[i,4] <- sum( F == 2 & M == 1 & O == 2, na.rm = TRUE)
      major[i,4] <- sum( F == 2 & M == 1 & O == 1, na.rm = TRUE)

      minor[i,5] <- sum( F == 1 & M == 1 & O == 1, na.rm = TRUE)
      major[i,5] <- minor[i,5]

      minor[i,6] <- 2*sum( F == 1 & M == 1 & O == 2, na.rm = TRUE)
      major[i,6] <- 2*sum( F == 1 & M == 1 & O == 0, na.rm = TRUE)
 #   })
  }
#  return(list( minor = rowSums(minor, na.rm = TRUE), major = rowSums(major, na.rm = TRUE), mendel = rowSums(mendel, na.rm = TRUE)))
  return(sum(minor, na.rm = TRUE))
})

setMethod("TransCount", signature( object = "SNPTrioExperiment", region = "GRangesList"), function( object, region ){
  minor <- numeric(length(region))
  for( i in 1:length(region) ){
    gr <- region[[i]]
    minor[i] <- TransCount(object = object, region = gr )
  }
  return(minor)
})

setMethod("relist.sgy", signature( object = "GRanges" ), function(object){
  results <- list()
  for( i in 1:length(object) ){
    results <- c( results, list(object[i]) )
  }
  return(GRangesList(results))
})

setMethod("setdiff", signature( x = "GRangesList", y = "GRanges" ), function( x, y, ... ){
  results <- list()
  for( i in 1:length(x) ){
    results <- c( results, list(setdiff( y, x[[i]], ...)))
  }
  return(GRangesList(results))
})

setMethod("ScanTrio", signature(object="SNPTrioExperiment", window = "GRanges", block = "GRanges"), function(object, window, block){
  window.list <- relist.sgy(window)
  window.out <- setdiff( window.list, block )  
  inside <- TransCount(object, window.list)
  outside <- TransCount(object, window.out)

  ## Need count of observeable transmissions.
  mat <- cbind(inside, outside)
  rownames(mat) <- names(window)
  return(mat)
})
