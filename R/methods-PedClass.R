setMethod("initialize", signature(.Object="PedClass"),
         function(.Object, ...){
           .Object
           callNextMethod(.Object, ... )
         })

setMethod( "PedClass", signature(object = "DataFrame"), function(object) {
  ped.DF <- with( as(object, "data.frame"), {
    DataFrame( famid = factor(famid),
              id = factor(id),
              fid = factor(fid),
              mid = factor(mid),
              sex = factor(sex),
              dx = factor(dx) )
  })
  new("PedClass", ped.DF)
})

setMethod("trios",  signature(object="PedClass"), function(object) {
  trio.df <- subset( as(object, "data.frame" ), !is.na(fid) & !is.na(mid) & !is.na(id) )
  trio.df[,c("id","fid","mid")]
})

## setMethod("ids",  signature(object="PedClass"), function(object){
##   ped.df <- pedigree(object)
##   with(ped.df, unique(c( id, fid, mid)))
## })
