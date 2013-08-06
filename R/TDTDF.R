setClass("TDTDF", contains="DataFrame", prototype = DataFrame( stat=numeric(), pval = numeric()))

setMethod("initialize", "TDTDF", function(.Object,  ... ){
  .Object <- callNextMethod()
  .Object
})

setMethod("show", signature(object="TDTDF"), function(.Object){
    o <- order(.Object$pval, decreasing=FALSE)
    .Object2 <- .Object[o,]
    callNextMethod(.Object2)
})

# The TDTDF object is initialized with:
# foo <- new("TDTDF", DataFrame(stat=rnorm(1e3),pval=runif(1e3)))
