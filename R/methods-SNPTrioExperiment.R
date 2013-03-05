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
  #GenomicRanges::show(object)
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
