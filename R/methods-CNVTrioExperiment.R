setMethod("initialize", "CNVTrioExperiment", function(.Object, pedigree, ... ){
  .Object@pedigree <- pedigree
  .Object <- callNextMethod()
  .Object
})

setMethod("CNVTrioExperiment", signature("SummarizedExperiment", "PedClass"), function(se, pedigree){
  new("CNVTrioExperiment", pedigree, se)
})

setMethod("show", signature(object="CNVTrioExperiment"), function(object){
  callNextMethod()
  cat("pedigree(", nrow(pedigree(object)), "): famid id fid mid sex dx\n", sep = "")
  cat("complete trios(", nrow(completeTrios(object)), "):\n", sep = "")
})

setMethod("pedigree",  signature(object="CNVTrioExperiment"), function(object) object@pedigree )

setMethod("completeTrios",  signature(object="CNVTrioExperiment"), function(object){
  ped <- pedigree(object)
  trios <- trios(ped)
  ids <- colnames(object)
  index <- (trios$id %in% ids) & (trios$fid %in% ids ) & (trios$mid %in% ids)
  return(data.frame( id = trios$id[index], mid = trios$mid[index], fid = trios$fid[index], stringsAsFactors = FALSE) )
})
