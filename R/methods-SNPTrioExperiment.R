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
  holger <- with( object, matrix(c(t(cbind( as(O,"matrix"), as(F,"matrix"), as(M,"matrix")))), byrow = TRUE, nrow = 3*nrow(O), ncol = ncol(O)))
  x <- ifelse( holger == 01, 0L, ifelse( holger == 02, 1L, ifelse( holger == 03, 2L, NA )))
  x
})

setMethod("TransCount", signature( object = "SNPTrioExperiment", gr = "GRanges"), function( object, gr ){
  ste <- object[subjectHits(findOverlaps(gr, rowData(object))),]
  show(ste)
  #ste <- object # temporary patch
  gtrio <- GenoTrio(ste)
  t <- numeric(6)
  t[1] <-   with( gtrio, sum( F == 1 & M == 2 & O == 2, na.rm = TRUE) )
  t[2] <-   with( gtrio, sum( F == 2 & M == 1 & O == 2, na.rm = TRUE) )
  t[3] <-   with( gtrio, sum( F == 2 & M == 3 & O == 3, na.rm = TRUE) )
  t[4] <-   with( gtrio, sum( F == 3 & M == 2 & O == 3, na.rm = TRUE) )
  t[5] <-   with( gtrio, sum( F == 2 & M == 2 & O == 2, na.rm = TRUE) )
  t[6] <- 2*with( gtrio, sum( F == 2 & M == 2 & O == 3, na.rm = TRUE) )
  return(sum(t))
})

setMethod("[", "SNPTrioExperiment", function( x, i, j, ..., drop = TRUE ){
  se <- as( x, "SummarizedExperiment" )
  if( !missing(j) ){
    se <- se[i,j]
  }else{
    se <- se[i,]
  }
  return(SNPTrioExperiment(se, pedigree(x)))
    
})
