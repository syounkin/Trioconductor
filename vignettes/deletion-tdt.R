### R code from vignette source 'deletion-tdt.Rnw'

###################################################
### code chunk number 1: options
###################################################
  options(width=75, continue = " ")
  library("Gviz")
  library("trioClasses")
  data("cnv", package = "trioClasses")
  data("pedigrees", package="CleftCNVAssoc")
  se <- SummarizedExperiment(assays = SimpleList(cnv = t(cnv.obj$cnv.mat)), colData = DataFrame(id=rownames(cnv.obj$cnv.mat), row.names = rownames(cnv.obj$cnv.mat)), rowData = cnv.obj$cmp.gr )
  beaty.trios <- MinimumDistance:::trios(beaty.pedigree)
  beaty.ped <- DataFrame(famid = do.call("rbind",strsplit(beaty.trios$O, "_" ))[,1], id = beaty.trios$O, fid = beaty.trios$F, mid = beaty.trios$M, sex = NA, dx = NA)
  ped <- PedClass(beaty.ped)
  fe <- FamilyExperiment(se, pedigree = ped )


###################################################
### code chunk number 2: FamilyExperiment
###################################################
  fe.parents <- fe[,colnames(fe)%in%parents(fe)]


###################################################
### code chunk number 3: freq-vec
###################################################
    freq.vec <- colSums(cnv(fe.parents))/nrow(cnv(fe.parents))


###################################################
### code chunk number 4: trioStates
###################################################
    trioAssay <- trioClasses:::TrioAssay(fe, type = "cnv")
    trioStates <- with(trioAssay, matrix( paste0(F,M,O), nrow = nrow(O), ncol = ncol(O)))
    dimnames(trioStates) <- dimnames(trioAssay$O)


###################################################
### code chunk number 5: table-list
###################################################
    table.list <- apply(trioStates, 2, "table")


###################################################
### code chunk number 6: transtab
###################################################
    trans.vec <- as( lapply( table.list, trioClasses:::trans.tab ), "numeric")
    head(table.list[which(trans.vec <= 0.05/length(trans.vec))])


