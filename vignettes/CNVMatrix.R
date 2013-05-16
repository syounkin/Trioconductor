### R code from vignette source './CNVMatrix.Rnw'

###################################################
### code chunk number 1: CNVMatrix.Rnw:3-4
###################################################
  options(width=70, continue = " ")


###################################################
### code chunk number 2: package
###################################################
library("trioClasses")
library("CleftCNVAssoc")
source("~/jhsph/R/packages/CleftCNVAssoc/vignettes/curated/make-data.R")

#compound.pitt <- grep("-", values(gr.pitt)$triostate)
#compound.beaty <- grep("-", values(gr.beaty)$triostate)

#gr.beaty <- gr.beaty[-compound.beaty]
#gr.pitt <- gr.pitt[-compound.pitt]

num.snp <- 10

gr.deletion.beaty <- gr.beaty[values(gr.beaty)$numsnp >= num.snp & values(gr.beaty)$cn %in% 0:1 ]
homos.beaty <- with(values(gr.deletion.beaty),cn==0)
#gr.deletion.beaty <- c(gr.deletion.beaty,gr.deletion.beaty[homos.beaty])
gr.deletion.pitt <- gr.pitt[values(gr.pitt)$numsnp >= num.snp & values(gr.pitt)$cn %in% 0:1 ]
homos.pitt <- with(values(gr.deletion.pitt),cn==0)
#gr.deletion.pitt <- c(gr.deletion.pitt,gr.deletion.pitt[homos.pitt])
gr.deletion.beaty.clean <- c( unlist(reduce(split(gr.deletion.beaty, values(gr.deletion.beaty)$id))), unlist(reduce(split(gr.deletion.beaty[homos.beaty], values(gr.deletion.beaty[homos.beaty])$id))))
gr.deletion.pitt.clean <- c( unlist(reduce(split(gr.deletion.pitt, values(gr.deletion.pitt)$id))), unlist(reduce(split(gr.deletion.pitt[homos.pitt], values(gr.deletion.pitt[homos.pitt])$id))))



###################################################
### code chunk number 3: cnvmatrix-beaty
###################################################
system.time( cnv.obj.beaty <- CNVMatrix( gr.deletion.beaty.clean, gr = gr.deletion.pitt.clean ) )


###################################################
### code chunk number 4: cnvmatrix-pitt
###################################################
system.time( cnv.obj.pitt <- CNVMatrix( gr.deletion.pitt.clean, gr = gr.deletion.beaty.clean ) )


###################################################
### code chunk number 5: savecnv
###################################################
save( cnv.obj.beaty, cnv.obj.pitt, file = "./../data/cnv.RData" )


