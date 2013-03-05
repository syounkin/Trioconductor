setMethod("initialize", "SNPTrioExperiment",
          function(.Object, pedigree, ... ){
            
            .Object@pedigree <- pedigree
            .Object <- callNextMethod()
             
            .Object
          })

#setMethod("SNPTrioExperiment", signature("SummarizedExperiment", "data.frame"), function(se, pedigree){
#  new("SNPTrioExperiment", se, pedigree)
#})

setMethod("show", signature(object="SNPTrioExperiment"), function(object){
  callNextMethod()
  #GenomicRanges::show(object)
  cat("pedigree(", nrow(pedigree(object)), "): FAMID ID FID MID SEX DX\n", sep = "")
  })

setMethod("geno", signature(object="SNPTrioExperiment"), function(object) assays(object)$geno )
setMethod("logR", signature(object="SNPTrioExperiment"), function(object) assays(object)$logR )
setMethod("baf",  signature(object="SNPTrioExperiment"), function(object) assays(object)$baf )


setMethod("pedigree",  signature(object="SNPTrioExperiment"), function(object) object@pedigree )

setMethod("completeTrios",  signature(object="SNPTrioExperiment"), function(object){
  ped.df <- pedigree(object)
  index.vec <- with( ped.df, id %in% colnames(snpEx) & fid %in% colnames(snpEx) & mid %in% colnames(snpEx) )
  ped.df[index.vec,]
})
