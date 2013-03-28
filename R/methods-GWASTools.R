getIndices <- function( gr,  snpAnnot ){
  # Note that this function assumes that snpAnnot is ordered by snpID!!!
  gr.gwas <- GRanges( seqnames = paste0("chr",getChromosome(snpAnnot)), ranges = IRanges(start = getPosition(snpAnnot), width = 1), strand = "*")
  snpIDs <- subjectHits(findOverlaps(gr, gr.gwas))
#  scanIDs <- getScanID(scanAnnot[with(as(scanAnnot,"data.frame"), scanName %in%  sub.vec)])
  return( snpIDs )
}
