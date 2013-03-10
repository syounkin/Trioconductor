setMethod("initialize", "SNPTrioExperiment",
          function(.Object, pedigree, ... ){
            
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

setMethod("TransCount", signature( object = "SNPTrioExperiment", gr = "GRanges"), function( object, gr ){
  minor <- matrix(0,nrow=length(gr), ncol = 6)
  major <- matrix(0,nrow=length(gr), ncol = 6)
  mendel <- matrix(0,nrow=length(gr), ncol = 4)
  for( i in 1:length(gr) ){
    gr.row <- gr[i]
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
  return(list( minor = rowSums(minor, na.rm = TRUE), major = rowSums(major, na.rm = TRUE), mendel = rowSums(mendel, na.rm = TRUE)))
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

