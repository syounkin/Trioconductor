setMethod("initialize", "SNPTrioExperiment",
          function(.Object,
                   assays = NA,
                   colData = NA,
                   rowData = NA,
                   pedigree = NA,
                   ...){
            .Object@pedigree <- pedigree
            callNextMethod(.Object, assays=assays, colData = colData, rowData = rowData, ...)
          })

setMethod("initialize", signature(.Object="PedClass"),
         function(.Object,
                  pedigree=data.frame(fid=character(),
                    mid=character(),
                    id=character(), stringsAsFactors=FALSE),
                  ...){
           callNextMethod(.Object, pedigree = pedigree, ...)
         })


setMethod("show", signature(object="SNPTrioExperiment"), function(object){
  GenomicRanges::show(object)
  cat("pedigree has rows =", dim(pedigree(object))[1], "\n")
  })

setMethod("geno", signature(object="SNPTrioExperiment"), function(object) assays(object)$geno )
setMethod("logR", signature(object="SNPTrioExperiment"), function(object) assays(object)$logR )
setMethod("baf",  signature(object="SNPTrioExperiment"), function(object) assays(object)$baf )
setMethod("pedigree",  signature(object="PedClass"), function(object) object@pedigree )

setMethod("pedigree",  signature(object="SNPTrioExperiment"), function(object) object@pedigree )

setMethod("completeTrios",  signature(object="SNPTrioExperiment"), function(object){
  ped.df <- pedigree(pedigree(object))
  index.vec <- with( ped.df, id %in% colnames(snpEx) & fid %in% colnames(snpEx) & mid %in% colnames(snpEx) )
  ped.df[index.vec,]
})
