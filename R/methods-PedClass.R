setMethod("initialize", signature(.Object="PedClass"),
         function(.Object, ...){
           .Object
           callNextMethod(.Object, ... )
         })

setMethod( "PedClass", signature(object = "data.frame"), function(object) new("PedClass", object) )

## setMethod("pedigree",  signature(object="PedClass"), function(object) object@pedigree )

## setMethod("ids",  signature(object="PedClass"), function(object){
##   ped.df <- pedigree(object)
##   with(ped.df, unique(c( id, fid, mid)))
## })
