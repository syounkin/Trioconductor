setMethod("initialize", signature(.Object="PedClass"),
         function(.Object,
                  pedigree=data.frame(fid=character(),
                    mid=character(),
                    id=character(), stringsAsFactors=FALSE),
                  ...){
           callNextMethod(.Object, pedigree = pedigree, ...)
         })

setMethod("pedigree",  signature(object="PedClass"), function(object) object@pedigree )

setMethod("ids",  signature(object="PedClass"), function(object){
  ped.df <- pedigree(object)
  with(ped.df, unique(c( id, fid, mid)))
})
