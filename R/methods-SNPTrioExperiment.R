setMethod("initialize", "SNPTrioExperiment",
          function(.Object,
                   assays = NA,
                   colData = NA,
                   rowData = NA,
                   Pedigree = NA,
                   ...){
            .Object@Pedigree <- Pedigree
            callNextMethod(.Object, assays=assays, colData = colData, rowData = rowData, ...)
          })

setMethod("show", signature(object="SNPTrioExperiment"), function(object){
  GenomicRanges::show(object)
  cat("pedigree has rows =", dim(object@Pedigree)[1], "\n")
  })

setMethod("geno", signature(object="SNPTrioExperiment"), function(object) assays(object)$geno )
setMethod("logR", signature(object="SNPTrioExperiment"), function(object) assays(object)$logR )
setMethod("baf",  signature(object="SNPTrioExperiment"), function(object) assays(object)$baf )
setMethod("Pedigree",  signature(object="SNPTrioExperiment"), function(object) object@Pedigree )
setMethod("completeTrios",  signature(object="SNPTrioExperiment"), function(object){
  index.vec <- with( Pedigree(object), id %in% colnames(snpEx) & fid %in% colnames(snpEx) & mid %in% colnames(snpEx) )
  Pedigree(object)[index.vec,]
})
