setMethod("initialize", signature(.Object="PedClass"),
         function(.Object, ...){
           .Object
           callNextMethod(.Object, ... )
         })

setValidity("PedClass", function(object){
	if(!all(c("famid", "id", "fid", "mid", "sex", "dx") %in% colnames(object))){
		msg <- "PedClass must have columns 'famid', 'id', 'fid', 'mid', 'sex', and 'dx'"
	} else msg <- TRUE
	msg
})

##setMethod( "PedClass", signature(object = "DataFrame"), function(object) {
##  ped.DF <- with( as(object, "data.frame"), {
##    DataFrame( famid = factor(famid),
##              id = factor(id),
##              fid = factor(fid),
##              mid = factor(mid),
##              sex = factor(sex),
##              dx = factor(dx) )
##  })
##  new("PedClass", ped.DF)
##})

## trivial constructor...assumes object is DataFrame
PedDataFrame <- function(){
    DataFrame(famid = factor(),
              id = factor(),
              fid = factor(),
              mid = factor(),
              sex = factor(),
              dx = factor())
}

## constructor
PedClass <- function(object){
	if(missing(object)) object <- PedDataFrame()
	new("PedClass", object)
}

setMethod("trios",  signature(object="PedClass"), function(object) {
  trio.df <- subset( as(object, "data.frame" ), !is.na(fid) & !is.na(mid) & !is.na(id) )
#  id.char <- as.character(trio.df$id)#,
  trio.df.2 <- data.frame( id = as.character(trio.df$id),
                          fid = as.character(trio.df$fid),
                          mid = as.character(trio.df$mid), stringsAsFactors=FALSE )
#                          )
  return(trio.df.2)
})

setMethod("parents",  signature(object="PedClass"), function(object){
    with(as(object,"data.frame"), unique(c(as.character(fid), as.character(mid))))
})

setMethod("offspring",  signature(object="PedClass"), function(object){
    with(as(object,"data.frame"), unique(as.character(id)))
})

setMethod("allSubjects",  signature(object="PedClass"), function(object){
    with(as(object,"data.frame"), unique(c(as.character(id), as.character(fid), as.character(mid))))
})
