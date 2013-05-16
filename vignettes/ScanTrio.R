### R code from vignette source './ScanTrio.Rnw'

###################################################
### code chunk number 1: ScanTrio.Rnw:3-4
###################################################
  options(width=70, continue = " ")


###################################################
### code chunk number 2: package
###################################################
library("trioClasses")
data("cleft.ts.pedigree")
data("8q24-european-all.sm")


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
ped <- PedClass(cleft.ts.pedigree)


###################################################
### code chunk number 6: fe
###################################################
(fe <- FamilyExperiment(se, pedigree = ped ))


###################################################
### code chunk number 7: missingness
###################################################
na.markers <- which(colSums(is.na(geno(fe)))>0) # Identify markers that have at least one NA call
(fe <- fe[-na.markers])                        # Removing these markers removes all missing data


###################################################
### code chunk number 8: scantrio
###################################################
wd <- 50e3
window <- GRanges( seqnames = seqnames(fe), ranges = IRanges(start = start(rowData(fe)), width = wd ) )
(scan.trio <- ScanTrio( object = fe[MAF(fe)<=0.01], window = window, block = range(rowData(fe))))


###################################################
### code chunk number 9: make-files-for-cpp
###################################################
trioClasses:::make.files.for.cpp(fe, fileroot = "./8q24")


