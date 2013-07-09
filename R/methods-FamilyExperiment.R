setMethod("initialize", "FamilyExperiment", function(.Object, pedigree, ... ){
  .Object@pedigree <- pedigree
  .Object <- callNextMethod()
  .Object
})

setMethod("FamilyExperiment", signature("SummarizedExperiment", "PedClass"), function(se, pedigree){
  new("FamilyExperiment", pedigree, se)
})

setMethod("show", signature(object="FamilyExperiment"), function(object){
  callNextMethod()
  cat("pedigree(", nrow(pedigree(object)), "): famid id fid mid sex dx\n", sep = "")
  cat("complete trios(", nrow(completeTrios(object)), "):\n", sep = "")
})

setMethod("geno", signature(object="FamilyExperiment"), function(object) {
  geno <- t(assays(object)$geno)
  geno
})

setMethod("cnv", signature(object="FamilyExperiment"), function(object) {
  geno <- t(assays(object)$cnv)
  geno
})

setMethod("logR", signature(object="FamilyExperiment"), function(object) assays(object)$logR )

setMethod("baf",  signature(object="FamilyExperiment"), function(object) assays(object)$baf )

setMethod("pedigree",  signature(object="FamilyExperiment"), function(object) object@pedigree )

setMethod("completeTrios",  signature(object="FamilyExperiment"), function(object){
  ped <- pedigree(object)
  trios <- trios(ped)
  ids <- colnames(object)
  index <- (trios$id %in% ids) & (trios$fid %in% ids ) & (trios$mid %in% ids)
  return(data.frame( id = trios$id[index], mid = trios$mid[index], fid = trios$fid[index], stringsAsFactors = FALSE) )
})

setMethod("TrioAssay",  signature(object="FamilyExperiment"), function(object, type = "geno"){
  ct <- completeTrios(object)
  if( type == "geno" )
    assay.mat <- as(geno(object),"numeric")
  if( type == "cnv" )
    assay.mat <- cnv(object)
  return(list( O = assay.mat[ct$id,], F = assay.mat[ct$fid,], M = assay.mat[ct$mid,]))
})

setMethod("ctcbind", signature( object = "list"), function( object ){
  holger <- with( object, matrix(c(t(cbind( as(F,"numeric"), as(M,"numeric"), as(O,"numeric")))), byrow = TRUE, nrow = 3*nrow(O), ncol = ncol(O)))
  holger
})

setMethod("[", c("FamilyExperiment","ANY","ANY"), function( x, i, j, ..., drop = TRUE ){
  se <- as( x, "SummarizedExperiment" )
  if( !missing(i) & !missing(j)){
    se <- se[i,j]
  }else if( !missing(i) & missing(j)){
    se <- se[i,]
  }else if( missing(i) & !missing(j)){
    se <- se[,j]
  }else {
    se <- se [,]
  }
  return(FamilyExperiment(se, pedigree(x)))
})

setMethod("parents",  signature(object="FamilyExperiment"), function(object){
  with(as(pedigree(object),"data.frame"), unique(c(as.character(fid), as.character(mid))))
})

setMethod("MAF", signature(object="FamilyExperiment"), function(object){
  with(col.summary(geno(object[,which(colnames(object) %in% parents(object))])),MAF)
})

setMethod("aTDT", signature(object="FamilyExperiment"), function(object){
  geno <- as( object, "matrix" )
  aTDT.fn(geno)
})

setMethod("aTDT", signature(object="matrix"), function(object){
  geno <- object
  aTDT.fn(geno)
})

setAs( from = "FamilyExperiment", to = "matrix", function(from){
  gtrio <- TrioAssay(from)
  holger <- ctcbind(gtrio)
  return(holger)
})

setMethod("TrivialGRangesList", signature( object = "GRanges" ), function(object){
  results <- list()
  for( i in 1:length(object) ){
    results <- c( results, list(object[i]) )
  }
  return(GRangesList(results))
})

setMethod("TransCount", signature( object = "FamilyExperiment", region = "GRanges"), function( object, region ){
  minor <- numeric(6)
  major <- numeric(6)
  mendel <- numeric(8)
    ste <- object[subjectHits(findOverlaps(region, rowData(object))),]
    gtrio <- TrioAssay(ste, type = "geno")

      Fa <- with(gtrio,F)
      Ma <- with(gtrio,M)
      Off <- with(gtrio,O)
      
      mendel[1] <- sum( Fa == 0 & Ma == 0 & ( Off == 1 | Off == 2 ), na.rm = TRUE)
      mendel[2] <- sum( Fa == 2 & Ma == 2 & ( Off == 0 | Off == 1 ), na.rm = TRUE)
      mendel[3] <- sum( Fa == 0 & Ma == 2 & ( Off == 0 | Off == 2 ), na.rm = TRUE)
      mendel[4] <- sum( Fa == 2 & Ma == 0 & ( Off == 0 | Off == 2 ), na.rm = TRUE)
    
      minor[1] <-  sum( Fa == 0 & Ma == 1 & Off == 1, na.rm = TRUE)
      major[1] <-  sum( Fa == 0 & Ma == 1 & Off == 0, na.rm = TRUE)
      mendel[5] <- sum( Fa == 0 & Ma == 1 & Off == 2, na.rm = TRUE)

      minor[2] <-  sum( Fa == 1 & Ma == 0 & Off == 1, na.rm = TRUE)
      major[2] <-  sum( Fa == 1 & Ma == 0 & Off == 0, na.rm = TRUE)
      mendel[6] <- sum( Fa == 1 & Ma == 0 & Off == 2, na.rm = TRUE)

      minor[3] <-  sum( Fa == 1 & Ma == 2 & Off == 2, na.rm = TRUE)
      major[3] <-  sum( Fa == 1 & Ma == 2 & Off == 1, na.rm = TRUE)
      mendel[7] <- sum( Fa == 1 & Ma == 2 & Off == 0, na.rm = TRUE)
  
      minor[4] <-  sum( Fa == 2 & Ma == 1 & Off == 2, na.rm = TRUE)
      major[4] <-  sum( Fa == 2 & Ma == 1 & Off == 1, na.rm = TRUE)
      mendel[8] <- sum( Fa == 2 & Ma == 1 & Off == 0, na.rm = TRUE)

      minor[5] <- sum( Fa == 1 & Ma == 1 & Off == 1, na.rm = TRUE)
      major[5] <- minor[5]

      minor[6] <- 2*sum( Fa == 1 & Ma == 1 & Off == 2, na.rm = TRUE)
      major[6] <- 2*sum( Fa == 1 & Ma == 1 & Off == 0, na.rm = TRUE)

  return(list( minor = sum(minor, na.rm = TRUE), major = sum(major, na.rm = TRUE), mendel = sum(mendel, na.rm = TRUE)))
})

setMethod("TransCount", signature( object = "FamilyExperiment", region = "GRangesList"), function( object, region ){
  minor <- numeric(length(region))
  major <- numeric(length(region))
  mendel <- numeric(length(region))  
  for( i in 1:length(region) ){
    gr <- region[[i]]
    trans <- TransCount(object = object, region = gr )
    minor[i] <- trans$minor
    major[i] <- trans$major
    mendel[i] <- trans$mendel
  }
  return(list( minor = minor, major = major, mendel = mendel))
})

setMethod("ScanTrio", signature(object="FamilyExperiment", window = "GRanges", block = "GRanges"), function(object, window, block){
  window.list <- TrivialGRangesList(window)
  trans.window <- TransCount(object, window.list)
  trans.block <- TransCount(object, block)
  df <-data.frame( minor.in = trans.window$minor, major.in = trans.window$major, mendel.in = trans.window$mendel, minor.out = trans.block$minor - trans.window$minor, major.out = trans.block$major - trans.window$major, mendel.out = trans.block$mendel - trans.window$mendel )
  rownames(df) <- names(window)
  with( df, {
    n <- (minor.in+minor.out+1)/(minor.in + minor.out + major.in + major.out + 2)
    y.in <- minor.in
    y.out <- minor.out
    n.in <- minor.in + major.in
    n.out <- minor.out + major.out    
    p.in <- (minor.in+1)/(n.in +2)
    p.out <- (minor.out+1)/(n.out +2)
    ## lr <- (p.in/n)^y.in*(p.out/n)^y.out*((1-p.in)/(1-n))^(n.in-y.in)*((1-p.out)/(1-n))^(n.out-y.out)
    loglr <- ifelse( p.in > p.out, y.in*log(p.in/n) + y.out*log(p.out/n) + (n.in-y.in)*log((1-p.in)/(1-n)) + (n.out-y.out)*log((1-p.out)/(1-n)), 0)
    lr <- exp(loglr)
    meta <- values(window)
    meta2 <- DataFrame( meta, lr = lr, minor.in = as.integer(minor.in), major.in = as.integer(major.in), minor.out = as.integer(minor.out), major.out = as.integer(major.out), mendel.in = as.integer(mendel.in), mendel.out = as.integer(mendel.out) )
    values(window) <- meta2
    return(window)
    ## DataFrame(lr = lr, minor.in = as.integer(minor.in), major.in = as.integer(major.in), minor.out = as.integer(minor.out), major.out = as.integer(major.out), mendel.in = as.integer(mendel.in), mendel.out = as.integer(mendel.out))
  })
})
