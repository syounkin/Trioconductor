setMethod("CNVMatrix",  signature(object="GRanges"), function(object){
  id.vec <- as.character(unique(values(object)$id))
  cmp <- disjoin(object)
  cnv <- matrix(NA, nrow = length(id.vec), ncol = length(cmp) )
  for ( i in 1:length(id.vec) ){
    cnv[i,] <- countOverlaps( cmp, object[values(object)$id == id.vec[i]] )
  }
  rownames(cnv) <- id.vec
  colnames(cnv) <- paste0("comp",1:length(cmp))
  names(cmp) <- colnames(cnv)
  return(list(cnv.mat = cnv, cmp.gr = cmp, gr = object ))
})


setMethod("intersect2", signature(object="GRanges"), function( object ){
  disjoin(object)[countOverlaps(disjoin(object),object) == length(object)]
})
