iTrio <- function(pedigreeData=Pedigree(), lrr=NULL, baf=NULL, geno=NULL){
 if(is.null(lrr) | is.null(baf) ){
    object <- new("iTrio", pedigree=pedigreeData)
 }else{
   data.mat <- lrr
   data.mat.2 <- baf
   father.names <- fatherNames(pedigreeData)
   mother.names <- motherNames(pedigreeData)
   offspring.names <- offspringNames(pedigreeData)
   nr <- nrow(data.mat)
   np <- length(offspring.names)
   lrrArray <- initializeBigArray("lrr", dim=c(nr, np, 3), vmode="integer")
   bafArray <- initializeBigArray("baf", dim=c(nr, np, 3), vmode="integer")   
   dimnames(lrrArray)[[1]] <- rownames(data.mat)
   dimnames(lrrArray)[[2]] <- offspring.names 
   dimnames(lrrArray)[[3]] <- c("F", "M", "O")
   dimnames(bafArray)[[1]] <- rownames(data.mat.2)
   dimnames(bafArray)[[2]] <- offspring.names 
   dimnames(bafArray)[[3]] <- c("F", "M", "O")
   lrrArray[,,"F"] <- data.mat[,father.names]
   lrrArray[,,"M"] <- data.mat[,mother.names]
   lrrArray[,,"O"] <- data.mat[,offspring.names]
   bafArray[,,"F"] <- data.mat.2[,father.names]
   bafArray[,,"M"] <- data.mat.2[,mother.names]
   bafArray[,,"O"] <- data.mat.2[,offspring.names]
   object <- new("iTrio", lrr=lrrArray, baf=bafArray,pedigree=pedigreeData)
 }
}

