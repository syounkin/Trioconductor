setMethod("initialize", "TrioSet",
          function(.Object,
                   assayData = assayDataNew(geno=geno),
                   #phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   #fatherPhenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   #motherPhenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   #annotation=character(),
                   #featureData=GenomeAnnotatedDataFrameFrom(assayData, annotation, genome=genome),
                   #experimentData=new("MIAME"),
                   #protocolData=phenoData[, integer(0)],
                   geno=array(NA, dim=c(0, 0, 3)),                   
                   #logRRatio=array(NA, dim=c(0, 0, 3)),
                   #BAF=array(NA, dim=c(0,0,3)),
                   pedigree=Pedigree(),
                   #mindist=NULL,
                   #genome=c("hg19", "hg18")
                    ...){
            .Object@pedigree <- pedigree
            #.Object@fatherPhenoData <- fatherPhenoData
            #.Object@motherPhenoData <- motherPhenoData
            callNextMethod(.Object,
                           assayData=assayData,
                           #phenoData=phenoData,
                           #fatherPhenoData=fatherPhenoData,
                           #motherPhenoData=motherPhenoData,
                           #featureData=featureData,
                           #experimentData=experimentData,
                           #annotation=annotation,
                           #protocolData=protocolData,
                           pedigree=pedigree,
                           #mindist=mindist,
                           #genome=match.arg(genome)
                            ...)
          })

setMethod("geno", "TrioSet",
          function(object) {
            assayDataElement(object, "geno")
          })
