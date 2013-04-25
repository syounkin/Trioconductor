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
  data("cnv", package = "trioClasses")
  cnv.beaty.obj <- cnv.obj
  data("cnv.pitt", package = "trioClasses")
  cnv.pitt.obj <- cnv.obj
  rm(cnv.obj)


###################################################
### code chunk number 2: FamilyExperiment
###################################################
  fe.beaty.parents <- fe.beaty[,colnames(fe.beaty)%in%parents(fe.beaty)]
  fe.pitt.parents <- fe.pitt[,colnames(fe.pitt)%in%parents(fe.pitt)]


###################################################
### code chunk number 3: freq-vec
###################################################
    freq.vec <- colSums(cnv(fe.beaty.parents))/nrow(cnv(fe.beaty.parents))
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
    table.list.beaty[which(trans.vec <= 0.05/sum(!is.na(trans.vec)))][[1]]
    
    binom.list.beaty <- lapply( table.list.beaty, trioClasses:::trans.rate )
    binom.list.pitt <- lapply( table.list.pitt, trioClasses:::trans.rate )

    trans.rate <- trans.vec #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    trans.rate.pitt <- trans.vec.pitt #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


###################################################
### code chunk number 7: gviz
###################################################
chr <- 8
gr.cnp <- cnv.beaty.obj$cmp.gr
gr.cnp.chr <- gr.cnp[as.logical(seqnames(gr.cnp) == paste0("chr", chr))]
    
p.vec <- trans.vec
p.vec.chr <- p.vec[as.logical(seqnames(gr.cnp) == paste0("chr", chr))]
freq.vec.chr <- freq.vec[as.logical(seqnames(gr.cnp) == paste0("chr", chr))]
trans.rate.chr <- trans.rate[as.logical(seqnames(gr.cnp) == paste0("chr", chr))]

TranscriptDb <- TxDb.Hsapiens.UCSC.hg18.knownGene
atrack <- AnnotationTrack(reduce(gr.cnp.chr), name = "CNP comp.", fill = "darkgreen")
gtrack <- GenomeAxisTrack()
dtrack <- DataTrack( range = gr.cnp.chr, data = -log10(p.vec.chr), type = "p", cex = 1, name = "-log10(p)", ylim = c(0,10), col = "red" )
dtrack2 <- DataTrack( range = gr.cnp.chr, data = ifelse(freq.vec.chr > 0.01, freq.vec.chr, NA ), type = "s", col = "blue", cex = 1, name = "Frequency", ylim = c(0,1), lwd = 2 )
dtrack3 <- DataTrack( range = gr.cnp.chr, data = trans.rate.chr, type = "s", col = "blue", cex = 1, name = "Transmission Rate", ylim = c(0,1), lwd = 2 )
itrack <- IdeogramTrack(genome = "hg18", chromosome = paste0("chr", chr), lty = 1, lwd = 1 )
grtrack <- GeneRegionTrack(TranscriptDb, genome="hg18", chromosome=chr, name="UCSC Known Transcripts", showId = TRUE)

fud <- 5e5/2
plotTracks(list(dtrack, dtrack2, dtrack3, grtrack, gtrack, itrack ),  background.panel = '#FFFEDB', background.title = 'darkblue', from = 39356825 - fud, to = 39497557 + fud)


###################################################
### code chunk number 8: gviz2
###################################################
chr <- 8
gr.cnp <- cnv.pitt.obj$cmp.gr
gr.cnp.chr <- gr.cnp[as.logical(seqnames(gr.cnp) == paste0("chr", chr))]
    
p.vec <- trans.vec.pitt
p.vec.chr <- p.vec[as.logical(seqnames(gr.cnp) == paste0("chr", chr))]
freq.vec.chr <- freq.pitt.vec[as.logical(seqnames(gr.cnp) == paste0("chr", chr))]
trans.rate.pitt.chr <- trans.rate.pitt[as.logical(seqnames(gr.cnp) == paste0("chr", chr))]

TranscriptDb <- TxDb.Hsapiens.UCSC.hg18.knownGene
atrack <- AnnotationTrack(reduce(gr.cnp.chr), name = "CNP comp.", fill = "darkgreen")
gtrack <- GenomeAxisTrack()
dtrack <- DataTrack( range = gr.cnp.chr, data = -log10(p.vec.chr), type = "p", cex = 1, name = "-log10(p)", ylim = c(0,10), col = "red" )
dtrack2 <- DataTrack( range = gr.cnp.chr, data = ifelse(freq.vec.chr > 0.01, freq.vec.chr, NA ), type = "s", col = "blue", cex = 1, name = "Frequency", ylim = c(0,1), lwd = 2 )
dtrack3 <- DataTrack( range = gr.cnp.chr, data = trans.rate.pitt.chr, type = "s", col = "blue", cex = 1, name = "Transmission Rate", ylim = c(0,1), lwd = 2 )
itrack <- IdeogramTrack(genome = "hg18", chromosome = paste0("chr", chr), lty = 1, lwd = 1 )
grtrack <- GeneRegionTrack(TranscriptDb, genome="hg18", chromosome=chr, name="UCSC Known Transcripts", showId = TRUE)

fud <- 5e5/2
plotTracks(list(dtrack, dtrack2, dtrack3, grtrack, gtrack, itrack ),  background.panel = '#FFFEDB', background.title = 'darkblue', from = 39356825 - fud, to = 39497557 + fud)


###################################################
### code chunk number 9: n-tests
###################################################
    n.tests <- sum(!is.na(trans.vec))


###################################################
### code chunk number 10: reduce-locus
###################################################
    (locus <- cnv.beaty.obj$cmp.gr[which(trans.vec <= 0.05/n.tests)])
    wd <- width(reduce(locus))/1e3


###################################################
### code chunk number 11: trans.tab
###################################################
    trioClasses:::trans.tab


###################################################
### code chunk number 12: pitt
###################################################
gr.deletion.pitt <- gr.pitt[values(gr.pitt)$numsnp >= 10 & values(gr.pitt)$cn %in% 0:1 ]
sum(countOverlaps( gr.deletion.pitt, reduce(locus) ))


###################################################
### code chunk number 13: hist
###################################################
par(bg = '#FFFEDB')
layout(matrix(1:2, nrow = 1, ncol = 2))
hist(trans.rate, breaks = 20, xlim = c(0,1), ylim = c(0,200), main = "Cleft Trios", xlab = "Transmission Rate", col='darkblue')
hist(trans.rate.pitt, breaks = 20, xlim = c(0,1), ylim = c(0,200), main = "Control Trios", xlab = "Transmission Rate", col='darkblue')


###################################################
### code chunk number 14: transmission-rate
###################################################
mean(trans.rate[subjectHits(findOverlaps(reduce(locus), cnv.beaty.obj$cmp.gr))])
mean(trans.rate.pitt[subjectHits(findOverlaps(reduce(locus), cnv.pitt.obj$cmp.gr))])


###################################################
### code chunk number 15: tdt-p
###################################################
mean(-log10(trans.vec[subjectHits(findOverlaps(reduce(locus), cnv.beaty.obj$cmp.gr))]))
mean(-log10(trans.vec.pitt[subjectHits(findOverlaps(reduce(locus), cnv.pitt.obj$cmp.gr))]))


###################################################
### code chunk number 16: frequencies
###################################################
freq.pitt.vec[subjectHits(findOverlaps(reduce(locus), cnv.pitt.obj$cmp.gr))]
freq.vec[subjectHits(findOverlaps(reduce(locus), cnv.beaty.obj$cmp.gr))]


###################################################
### code chunk number 17: ci-beaty
###################################################
est.list.beaty <- lapply(binom.list.beaty, FUN = trioClasses:::get.est)
ci.list.beaty <- lapply(binom.list.beaty, FUN = trioClasses:::get.ci)
testable.beaty <- !is.na(ci.list.beaty)

est.vec.beaty <- as(est.list.beaty, "numeric")[testable.beaty]
cnv.beaty.obj$cmp.gr[testable.beaty]


###################################################
### code chunk number 18: ci-pitt
###################################################
est.list.pitt <- lapply(binom.list.pitt, FUN = trioClasses:::get.est)
ci.list.pitt <- lapply(binom.list.pitt, FUN = trioClasses:::get.ci)
testable.pitt <- !is.na(ci.list.pitt)
ci.mat.pitt <- matrix(unlist(ci.list.pitt[testable.pitt]), nrow = sum(testable.pitt), ncol = 2, byrow = TRUE)
est.vec.pitt <- as(est.list.pitt, "numeric")[testable.pitt]
cnv.pitt.obj$cmp.gr[testable.pitt]


###################################################
### code chunk number 19: trans-ci
###################################################
par(bg = '#FFFEDB')
plot(1, type = "n", ylim = c(0,length(est.vec.pitt)), xlim = c(0,1), xlab = "Transmission Rate", ylab = "", main = "Cleft Trios Transmission Rate 95 percent CI\n40 0/1 mating pairs required" )
polygon( y = c(1:length(est.vec.beaty), length(est.vec.beaty):1 ), x = c(ci.mat.beaty[,1],rev(ci.mat.beaty[,2])), col = "darkblue", border = NA)
lines( x = rep(0.5,2), y = c(0,length(est.vec.beaty)), lty = 2, lwd = 2)
#points( x = est.vec, y = (1:length(est.vec))-1, pch = 20, cex = 0.5 )
par(bg = '#FFFEDB')
plot(1, type = "n", ylim = c(0,length(est.vec.pitt)), xlim = c(0,1), xlab = "Transmission Rate", ylab = "", main = "Control Trios Transmission Rate 95 percent CI\n40 0/1 mating pairs required" )
polygon( y = c(1:length(est.vec.pitt), length(est.vec.pitt):1 ), x = c(ci.mat.pitt[,1],rev(ci.mat.pitt[,2])), col = "darkblue", border = NA)
lines( x = rep(0.5,2), y = c(0,length(est.vec.pitt)), lty = 2, lwd = 2)
#points( x = est.vec, y = (1:length(est.vec))-1, pch = 20, cex = 0.5 )


