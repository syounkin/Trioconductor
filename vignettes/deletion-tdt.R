### R code from vignette source 'deletion-tdt.Rnw'

###################################################
### code chunk number 1: options
###################################################
  options(width=75, continue = " ")
  library("Gviz")
  library("trioClasses")
  library("TxDb.Hsapiens.UCSC.hg18.knownGene")
  data("fe", package = "trioClasses")
  data("pedigrees", package="CleftCNVAssoc")
  data("penncnvjoint", package = "CleftCNVAssoc")
#  data("cnv", package = "trioClasses")
#  cnv.obj.beaty <- cnv.obj
#  data("cnv.pitt", package = "trioClasses")
#  cnv.pitt.obj <- cnv.obj
#  rm(cnv.obj)


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
### code chunk number 6: transtab
###################################################
    trans.vec <- as( lapply( table.list.beaty, trioClasses:::trans.tab ), "numeric")
    trans.vec.pitt <- as( lapply( table.list.pitt, trioClasses:::trans.tab ), "numeric")
#    table.list.beaty[which(trans.vec <= 0.05/sum(!is.na(trans.vec)))][[1]]
    
    binom.list.beaty <- lapply( table.list.beaty, trioClasses:::trans.rate )
    binom.list.pitt <- lapply( table.list.pitt, trioClasses:::trans.rate )

#    trans.rate <- trans.vec #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#    trans.rate.pitt <- trans.vec.pitt #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


###################################################
### code chunk number 7: ci-beaty
###################################################
est.list.beaty <- lapply(binom.list.beaty, FUN = trioClasses:::get.est)
ci.list.beaty <- lapply(binom.list.beaty, FUN = trioClasses:::get.ci)
testable.beaty <- !is.na(ci.list.beaty)
ci.mat.beaty <- matrix(unlist(ci.list.beaty[testable.beaty]), nrow = sum(testable.beaty), ncol = 2, byrow = TRUE)
est.vec.beaty <- as(est.list.beaty, "numeric")[testable.beaty]
gr.cnp.beaty <- cnv.obj.beaty$cmp.gr[testable.beaty]
freq.beaty.vec <- freq.beaty.vec[testable.beaty]


###################################################
### code chunk number 8: ci-pitt
###################################################
est.list.pitt <- lapply(binom.list.pitt, FUN = trioClasses:::get.est)
ci.list.pitt <- lapply(binom.list.pitt, FUN = trioClasses:::get.ci)
testable.pitt <- !is.na(ci.list.pitt)
ci.mat.pitt <- matrix(unlist(ci.list.pitt[testable.pitt]), nrow = sum(testable.pitt), ncol = 2, byrow = TRUE)
est.vec.pitt <- as(est.list.pitt, "numeric")[testable.pitt]
gr.cnp.pitt <- cnv.obj.pitt$cmp.gr[testable.pitt]
freq.pitt.vec <- freq.pitt.vec[testable.pitt]


###################################################
### code chunk number 9: trans-ci
###################################################
par(bg = '#FFFEDB')
plot(1, type = "n", ylim = c(0,max(length(est.vec.beaty),length(est.vec.pitt))), xlim = c(0,2), xlab = "Transmission Rate", ylab = "", main = "Transmission Rate 95 percent CI\n60 0/1 mating pairs required", axes = FALSE )
polygon( y = c(1:length(est.vec.beaty), length(est.vec.beaty):1 ), x = c(ci.mat.beaty[,1],rev(ci.mat.beaty[,2])), col = "darkblue", border = NA)

polygon( y = c(1:length(est.vec.pitt), length(est.vec.pitt):1 ), x = c(ci.mat.pitt[,1],rev(ci.mat.pitt[,2]))+1, col = "orange", border = NA)
#lines( x = rep(0.5,2), y = c(0,length(est.vec.beaty)), lty = 2, lwd = 2)
lines( x = rep(0.5,2), y = c(0,max(length(est.vec.beaty),length(est.vec.pitt))), lty = 2, lwd = 2)
lines( x = rep(0.5,2)+1, y = c(0,max(length(est.vec.beaty),length(est.vec.pitt))), lty = 2, lwd = 2)
lines( x = rep(1,2), y = c(0,max(length(est.vec.beaty),length(est.vec.pitt))), lty = 3, lwd = 2)
axis(1, at = c(0.5,1.5), labels = c(0.5,0.5) )
axis(2)


###################################################
### code chunk number 10: n-tests
###################################################
    n.tests <- sum(!is.na(trans.vec))


###################################################
### code chunk number 11: reduce-locus
###################################################
    (locus <- cnv.obj.beaty$cmp.gr[which(trans.vec <= 0.05/n.tests)])
    wd <- width(reduce(locus))/1e3


###################################################
### code chunk number 12: trans.tab
###################################################
    trioClasses:::trans.tab


###################################################
### code chunk number 13: pitt
###################################################
gr.deletion.pitt <- gr.pitt[values(gr.pitt)$numsnp >= 10 & values(gr.pitt)$cn %in% 0:1 ]
sum(countOverlaps( gr.deletion.pitt, reduce(locus) ))


###################################################
### code chunk number 14: hist
###################################################
par(bg = '#FFFEDB')
layout(matrix(1:2, nrow = 1, ncol = 2))
hist(est.vec.beaty, breaks = 20, xlim = c(0,1), ylim = c(0,200), main = "Cleft Trios", xlab = "Transmission Rate", col='darkblue')
hist(est.vec.pitt, breaks = 20, xlim = c(0,1), ylim = c(0,200), main = "Control Trios", xlab = "Transmission Rate", col='darkblue')


###################################################
### code chunk number 15: gviz
###################################################
chr <- 8
#gr.cnp.beaty <- cnv.obj.beaty$cmp.gr
gr.cnp.beaty.chr <- gr.cnp.beaty[as.logical(seqnames(gr.cnp.beaty) == paste0("chr", chr))]

#gr.cnp.pitt <- cnv.obj.pitt$cmp.gr
gr.cnp.pitt.chr <- gr.cnp.pitt[as.logical(seqnames(gr.cnp.pitt) == paste0("chr", chr))]

#p.vec <- trans.vec
est.vec.beaty.chr <- est.vec.beaty[as.logical(seqnames(gr.cnp.beaty) == paste0("chr", chr))]
freq.vec.beaty.chr <- freq.beaty.vec[as.logical(seqnames(gr.cnp.beaty) == paste0("chr", chr))]
est.vec.pitt.chr <- est.vec.pitt[as.logical(seqnames(gr.cnp.pitt) == paste0("chr", chr))]
freq.vec.pitt.chr <- freq.pitt.vec[as.logical(seqnames(gr.cnp.pitt) == paste0("chr", chr))]
#trans.rate.chr <- trans.rate[as.logical(seqnames(gr.cnp) == paste0("chr", chr))]

TranscriptDb <- TxDb.Hsapiens.UCSC.hg18.knownGene
atrack <- AnnotationTrack(, name = "CNP comp.", fill = "darkgreen")
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg18", chromosome = paste0("chr", chr), lty = 1, lwd = 1 )
grtrack <- GeneRegionTrack(TranscriptDb, genome="hg18", chromosome=chr, name="UCSC Known Genes", geneSymbols = TRUE, collapseTranscripts=TRUE, showId = TRUE)

dtrack <- DataTrack( range = gr.cnp.beaty.chr, data = rbind(est.vec.beaty.chr, freq.vec.beaty.chr), type = "b", cex = 1, name = "Cleft Trios", ylim = c(0,1), col = c("red","darkblue"), pch = 20 )
dtrack2 <- DataTrack( range = gr.cnp.pitt.chr, data = rbind(est.vec.pitt.chr, freq.vec.pitt.chr), type = "b", cex = 1, name = "Control Trios", ylim = c(0,1), col = c("red","darkblue"), pch = 20 )
fud <- 2e5/2
plotTracks( list(dtrack, dtrack2, grtrack, gtrack, itrack ),  groups = c("trans","freq"), background.panel = '#FFFEDB', background.title = 'darkblue', from = 39356825 - fud, to = 39497557 + fud)

#dtrack2 <- DataTrack( range = gr.cnp.pitt.chr, data = est.vec.pitt.chr, type = "b", cex = 1, name = "Trans (Controls)", ylim = c(0,1), col = "darkblue", pch = 20 )
#dtrack3 <- DataTrack( range = gr.cnp.chr, data = ifelse(freq.vec.chr > 0.01, freq.vec.chr, NA ), type = "s", col = "blue", cex = 1, name = "Frequency", ylim = c(0,1), lwd = 2 )
#dtrack3 <- DataTrack( range = gr.cnp.chr, data = trans.rate.chr, type = "s", col = "blue", cex = 1, name = "Transmission Rate", ylim = c(0,1), lwd = 2 )
#itrack <- IdeogramTrack(genome = "hg18", chromosome = paste0("chr", chr), lty = 1, lwd = 1 )
#grtrack <- GeneRegionTrack(TranscriptDb, genome="hg18", chromosome=chr, name="UCSC Known Transcripts", showId = TRUE)

#fud <- 5e5/2
#plotTracks(list(dtrack,  grtrack, gtrack, itrack ),  background.panel = '#FFFEDB', background.title = 'darkblue', from = 39356825 - fud, to = 39497557 + fud)


###################################################
### code chunk number 16: transmission-rate
###################################################
mean(est.vec.beaty[subjectHits(findOverlaps(reduce(locus), cnv.obj.beaty$cmp.gr))])
mean(est.vec.pitt[subjectHits(findOverlaps(reduce(locus), cnv.obj.pitt$cmp.gr))])


###################################################
### code chunk number 17: tdt-p
###################################################
mean(-log10(trans.vec[subjectHits(findOverlaps(reduce(locus), cnv.obj.beaty$cmp.gr))]))
mean(-log10(trans.vec.pitt[subjectHits(findOverlaps(reduce(locus), cnv.obj.pitt$cmp.gr))]))


###################################################
### code chunk number 18: frequencies
###################################################
freq.pitt.vec[subjectHits(findOverlaps(reduce(locus), cnv.obj.pitt$cmp.gr))]
freq.vec[subjectHits(findOverlaps(reduce(locus), cnv.obj.beaty$cmp.gr))]


