### R code from vignette source './inherited-deletions.Rnw'

###################################################
### code chunk number 1: options
###################################################
  options(width=75, continue = " ")
  library("Gviz")
  library("trioClasses")
  library("TxDb.Hsapiens.UCSC.hg18.knownGene")
##  data("fe-25", package = "trioClasses")
  data("fe", package = "trioClasses")
  data("pedigrees", package="CleftCNVAssoc")
  data("penncnvjoint", package = "CleftCNVAssoc")
  data("cnv", package = "trioClasses")
##data("cnv-25", package = "trioClasses")


###################################################
### code chunk number 2: FamilyExperiment
###################################################
  fe.beaty.parents <- fe.beaty[,colnames(fe.beaty)%in%parents(fe.beaty)]
  fe.pitt.parents <- fe.pitt[,colnames(fe.pitt)%in%parents(fe.pitt)]


###################################################
### code chunk number 3: freq-vec
###################################################
    freq.beaty.vec <- colSums(cnv(fe.beaty.parents))/nrow(cnv(fe.beaty.parents))
    freq.pitt.vec <- colSums(cnv(fe.pitt.parents))/nrow(cnv(fe.pitt.parents))


###################################################
### code chunk number 4: trioStates
###################################################
    trioAssay.beaty <- trioClasses:::TrioAssay(fe.beaty, type = "cnv")
    trioStates.beaty <- with(trioAssay.beaty, matrix( paste0(F,M,O), nrow = nrow(O), ncol = ncol(O)))
    dimnames(trioStates.beaty) <- dimnames(trioAssay.beaty$O)
    trioAssay.pitt <- trioClasses:::TrioAssay(fe.pitt, type = "cnv")
    trioStates.pitt <- with(trioAssay.pitt, matrix( paste0(F,M,O), nrow = nrow(O), ncol = ncol(O)))
    dimnames(trioStates.pitt) <- dimnames(trioAssay.pitt$O)


###################################################
### code chunk number 5: table-list
###################################################
    table.list.beaty <- apply(trioStates.beaty, 2, "table")
    table.list.pitt <- apply(trioStates.pitt, 2, "table")


###################################################
### code chunk number 6: TU
###################################################
TU.mat.beaty <- matrix(unlist(lapply(table.list.beaty, trioClasses:::CountTU)), nrow = length(table.list.beaty), ncol = 2, byrow = TRUE )
TU.mat.pitt <- matrix(unlist(lapply(table.list.pitt, trioClasses:::CountTU)), nrow = length(table.list.pitt), ncol = 2, byrow = TRUE )
TU.mat <- cbind(TU.mat.beaty, TU.mat.pitt)
testable <- which(   (rowSums(TU.mat[,1:2])>=25) & (rowSums(TU.mat[,3:4])>=25) )
TU.mat <- TU.mat[testable,]
rownames(TU.mat) <- names(table.list.beaty)[testable]
colnames(TU.mat) <- c("T.case","U.case","T.con","U.con")
DF <- DataFrame(TU.mat, rowData(fe.beaty)[testable])
colnames(DF) <- c(colnames(TU.mat),"grange")


###################################################
### code chunk number 7: hist
###################################################
hist(trans.vec <- rowSums(TU.mat[,c(1,3)])/rowSums(TU.mat), breaks = 20)


###################################################
### code chunk number 8: fish
###################################################
fish.list <- apply(TU.mat,1,trioClasses:::TU.fish)
p.vec <- unlist( lapply( fish.list, function(obj) return(obj$p.value) ) )
DF <- DataFrame(DF,p.vec)


###################################################
### code chunk number 9: length1
###################################################
c(length(DF$grange),length(reduce(DF$grange)))


###################################################
### code chunk number 10: qqplot
###################################################
n <- nrow(DF)
plot( -log10((1:n)/n), -log10(DF$p.vec[order(DF$p.vec)]), xlim = xlim <- c(0,5), ylim = xlim)
lines( c(0,xlim[2]), c(0,xlim[2]), lty = 3 )


###################################################
### code chunk number 11: badlocus
###################################################
badloci.gr <- reduce(rowData(fe.beaty[testable]))[which(countOverlaps(reduce(rowData(fe.beaty[testable])),rowData(fe.beaty[testable]))>10)]
bad <- subjectHits(findOverlaps( badloci.gr,DF$grange))
DF <- DF[-bad,]
TU.mat <- TU.mat[-bad,]


###################################################
### code chunk number 12: length2
###################################################
c(length(DF$grange),length(reduce(DF$grange)))


###################################################
### code chunk number 13: hist2
###################################################
hist(trans.vec <- rowSums(TU.mat[,c(1,3)])/rowSums(TU.mat), breaks = 10)


###################################################
### code chunk number 14: qqplot2
###################################################
n <- nrow(DF)
plot( -log10((1:n)/n), -log10(DF$p.vec[order(DF$p.vec)]), xlim = xlim <- c(0,5), ylim = xlim)
lines( c(0,xlim[2]), c(0,xlim[2]), lty = 3 )


