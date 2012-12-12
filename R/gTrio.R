gTrio <- function(pedigreeData=Pedigree(), lrr=NULL, baf=NULL, geno=NULL, map = data.frame()){
 if(is.null(geno)){
    object <- new("gTrio", pedigree=pedigreeData)
 }else{
   data.mat <- geno
   father.names <- fatherNames(pedigreeData)
   mother.names <- motherNames(pedigreeData)
   offspring.names <- offspringNames(pedigreeData)
   nr <- nrow(data.mat)
   np <- length(offspring.names)
   genoArray <- initializeBigArray("geno", dim=c(nr, np, 3), vmode="integer")
   dimnames(genoArray)[[1]] <- rownames(data.mat)
   dimnames(genoArray)[[2]] <- offspring.names 
   dimnames(genoArray)[[3]] <- c("F", "M", "O")
   genoArray[,,"F"] <- data.mat[,father.names]
   genoArray[,,"M"] <- data.mat[,mother.names]
   genoArray[,,"O"] <- data.mat[,offspring.names]
   object <- new("gTrio", geno=genoArray, pedigree=pedigreeData, map = map)
 }
}

##   object <- new("TrioSet",
##      #           BAF=bafArray,
##       #          logRRatio=logRArray,
##                 geno=genoArray,
## #                phenoData=phenoData,
##  #               fatherPhenoData=fatherPhenoData,
##   #              motherPhenoData=motherPhenoData,
##                 pedigree=pedigreeData
##                 #featureData=fD,
##    #             mindist=mindist,
##     #            genome=match.arg(genome))


 # np <- nrow(trios(pedigreeData))
 # trio.names <- array(NA, dim=c(length(offspringNames(pedigreeData)), 1, 3))
 # dimnames(trio.names) <- list(offspringNames(pedigreeData), "sampleNames", c("F", "M", "O"))
 # trio.names[, "sampleNames", ] <- as.matrix(trios(pedigreeData))
 #father.index <- match(father.names, colnames(data.mat))
 #mother.index <- match(mother.names,colnames(data.mat))
 #offspring.index <- match(offspring.names,colnames(data.mat))


                                        #  if(length(father.index) == 0) stop("father ids in pedigree do not match any of the column names of the data.mat matrix")
                                        #  if(length(mother.index) == 0) stop("mother ids in pedigree do not match any of the column names of the data.mat matrix")

 
#  if(length(offspring.index) == 0) stop("offspring ids in pedigree do not match any of the column names of the data.mat matrix")
 # bafArray <- initializeBigArray("baf", dim=c(nr, np, 3), vmode="integer")
 # logRArray <- initializeBigArray("lrr", dim=c(nr, np, 3), vmode="integer")
 # logRArray[,,"F"] <- lrr[, father.index]
 # logRArray[,,"M"] <- lrr[, mother.index]
 # logRArray[,,"O"] <- lrr[, offspring.index]
 # bafArray[,,"F"] <- baf[, father.index]
 # bafArray[,,"M"] <- baf[, mother.index]
 # bafArray[,,"O"] <- baf[, offspring.index]
  #genoArray[,,"M"] <- geno[, mother.index]
  #genoArray[,,"O"] <- geno[, offspring.index]
  ## if(nrow(pedigreeData) > 0){
  ##   if(!missing(sample.sheet)){
  ##     if(is.null(row.names)){
  ##       row.names <- rownames(sample.sheet)
  ##     }
  ##     if(!all(row.names %in% allNames(pedigreeData))){
  ##       stop("There are row.names for sample.sheet not in the pedigree object")
  ##     }
  ##     phenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE,
  ##                                         sample.sheet=sample.sheet,
  ##                                         which="offspring",
  ##                                         row.names=row.names)
  ##     fatherPhenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE,
  ##                                               sample.sheet=sample.sheet,
  ##                                               which="father",
  ##                                               row.names=row.names)
  ##     motherPhenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE,
  ##                                               sample.sheet=sample.sheet,
  ##                                               which="mother",
  ##                                               row.names=row.names)
  ##   } else {
      ## phenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE, which="offspring")
      ## fatherPhenoData <- annotatedDataFrameFrom(pedigreeData, FALSE, which="father")
      ## motherPhenoData <- annotatedDataFrameFrom(pedigreeData, FALSE, which="mother")
  ##   }
  ## }


  #if(missing(lrr) | missing(baf)){
  #  object <- new("TrioSet",
   #               pedigree=pedigreeData)
   # return(object)
  #} else{
   # if(ncol(lrr) > 0 & nrow(pedigreeData)==0)
    #  stop("pedigreeData has zero rows")
  #}
  ## if(!missing(lrr) & !missing(baf)){
  ##   if(!identical(rownames(lrr), rownames(baf)))
  ##     stop("rownames of lrr and baf are not identical")
  ##   if(!identical(dim(lrr), dim(baf)))
  ##     stop("lrr and baf must have the same dimension")
  ##   if(!(is(lrr[1,1], "integer") & is(baf[1,1], "integer"))){
  ##     stop("rr and baf must be integers. Use integerMatrix(x, scale=100) to transform log R ratios and integerMatrix(x, scale=1000) for B allele frequencies")
  ##   }
  ## }
  #if(missing(featureData)){
  #  if(missing(cdfname)) stop("If featureData is not supplied, a valid cdfname must be provided for feature annotation")
   # featureData <- GenomeAnnotatedDataFrameFrom(lrr, cdfname, genome=match.arg(genome))
    #fD <- featureData[order(chromosome(featureData), position(featureData)), ]
    #rm(featureData); gc()
  #} else {
   # if(!is(featureData, "AnnotatedDataFrame")) stop("featureData must be an AnnotatedDataFrame or a GenomeAnnotatedDataFrame")
    #fD <- featureData
  #}
#  is.present <- sampleNames(fD) %in% rownames(lrr)
 # if(!all(is.present)) fD <- fD[is.present, ]
 # if(!is.null(rownames(lrr))){
  #  index <- match(sampleNames(fD), rownames(lrr))
   # if(length(index) == 0) {
    #  if(!missing(cdfname)){
     #   msg <- paste("rownames for log R ratios do not match feature ids with annotation package ", cdfname)
      #  stop(msg)
 #     }
#}
#    lrr <- lrr[index, ]
#    baf <- baf[index, ]
#    stopifnot(all(identical(rownames(lrr), sampleNames(fD))))
#}
#  if(!drop){
 #   dimnames(bafArray)[c(1,2)] <- dimnames(logRArray)[c(1,2)] <- list(sampleNames(fD), colnames(lrr)[offspring.index])
# }
