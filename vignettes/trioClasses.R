### R code from vignette source './trioClasses.Rnw'

###################################################
### code chunk number 1: trioClasses.Rnw:3-4
###################################################
  options(width=70, continue = " ")


###################################################
### code chunk number 2: package
###################################################
#Rprof(filename = "trioClasses.Rprof")
library("trioClasses")
data("sample")

sge <- as(Sys.getenv("SGE_TASK_ID"),"integer")

cat("The SGE task ID is ", sge, "\n", sep = "")

if( sge%%3 == 1 ){
  data("8q24-european-all.sm")
}else if( sge%%3 == 2) {
  data("8q24-chinese-all.sm")
}else if( sge%%3 == 0 ){
  data("8q24-filipino-all.sm")
}else{
  cat("Problem with SGE task ID\n")
}

window.vec <- c(50,100,250)*1e3
k <- length(window.vec)

# hard-coded for 3 window lengths
if( sge %in% 1:k ){
  windowlength <- window.vec[1]
}else if( sge %in% (k+1):(2*k) ){
  windowlength <- window.vec[2]
}else if( sge %in% (2*k+1):(3*k)){
  windowlength <- window.vec[3]
}else{
  cat("Problem with SGE task ID\n")
}
  

###################################################
### code chunk number 3: prelim
###################################################
pos <- as(do.call("rbind", strsplit(colnames(sm),split = ":") )[,2], "integer" )
chr <- do.call("rbind", strsplit(colnames(sm),split = ":") )[,1]
gr <- GRanges( seqnames = chr, ranges = IRanges(start = pos, width = 1), strand = "*" )
names(gr) <- colnames(sm)
col.DF <- DataFrame( id = rownames(sm) )
rownames(col.DF) <- rownames(sm)


###################################################
### code chunk number 4: data
###################################################
se <- SummarizedExperiment(assays = SimpleList(geno = t(sm)), colData = col.DF, rowData = gr )


###################################################
### code chunk number 5: ped
###################################################
ped <- PedClass(ped.DF)


###################################################
### code chunk number 6: ste
###################################################
(ste <- FamilyExperiment(se, pedigree = ped ))


###################################################
### code chunk number 7: gtrio (eval = FALSE)
###################################################
## gtrio <- trioClasses:::TrioAssay(ste)
## index <- with(gtrio, !( is.na(F) | is.na(M) | is.na(O) ) )
## F <- as(with(gtrio,F), "numeric")[index]
## M <- as(with(gtrio,M), "numeric")[index]
## O <- as(with(gtrio,O), "numeric")[index]
## table(paste0(F,M,O))


###################################################
### code chunk number 8: sterare
###################################################
(ste.rare <- ste[!(MAF(ste) >= 0.01 | is.na(MAF(ste)))])


###################################################
### code chunk number 9: window-by-snv
###################################################
## n.win <- 100
## n.snp <- ceiling(length(rowData(ste.rare))/n.win)
## p <- length(rowData(ste.rare))
## i.vec <- seq(1,p,by = n.snp )
## j.vec <- seq(n.snp,p+n.snp,by = n.snp )
## j.vec[length(j.vec)] <- p
## #n <- 1000
start <- start(range(rowData(ste.rare)))
end <- end(range(rowData(ste.rare)))
## #wd <- ceiling((end-start)/n)
## #wd <- 100e3
## #start.vec <- seq(start,end+wd,by = wd)
## window <- GRanges( seqnames = levels(seqnames(rowData(ste.rare))), ranges = IRanges( start = start(rowData(ste.rare))[i.vec], end = start(rowData(ste.rare))[j.vec] ), strand = "*" )
## window <- window[-length(window)]


###################################################
### code chunk number 10: window (eval = FALSE)
###################################################
## n.snp
## table(countOverlaps( window, rowData(ste.rare)))


###################################################
### code chunk number 11: scantrio
###################################################
window <- rowData(ste.rare) + windowlength
system.time( scan.trio <- ScanTrio( object = ste.rare, window = window, block = range(rowData(ste.rare))))
scan.trio


###################################################
### code chunk number 12: plotscantrio
###################################################
pdf(file = paste0("./trioC-", sge, ".pdf"))
with( as( scan.trio, "data.frame" ), {
  plot( rowSums(cbind(start(window),end(window)))/2, log(lr), pch = 20, type = "b", lty = 3, axes = FALSE, xlab = "MB", ylab = "log(LR)")
  axis(2)
axis(1, at = at <- seq(start,end,length.out=4), labels = round(at/1e6,1) )
})
dev.off()

###################################################
### code chunk number 13: save
###################################################
save(scan.trio, file = paste0("./../data/scan-trio-",sge,".RData") )
