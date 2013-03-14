setMethod("CNVMatrix",  signature(object="GRanges"), function(object){
  id <- as.character(unique(values(object)$id))
  deletion <- object[values(object)$cn < 3 ]
  cmp <- disjoin(deletion)
  cnv <- matrix(NA, nrow = length(id), ncol = length(cmp) )
  for ( i in 1:length(id) ){
    cnv[i,] <- countOverlaps( cmp, deletion[with(values(deletion),id == id[i])])
  }
  rownames(cnv) <- id
  colnames(cnv) <- paste0("comp",1:length(cmp))
  names(cmp) <- colnames(cnv)
  return(list(cnv.mat = cnv, cmp.gr = cmp ))
})
