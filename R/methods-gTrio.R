setMethod("initialize", "gTrio",
          function(.Object,
                   assayData = assayDataNew(geno=geno),
                   geno=array(NA, dim=c(0, 0, 3)),                   
                   pedigree=Pedigree(), ...)
          {
            .Object@pedigree <- pedigree
            callNextMethod(.Object,
            assayData=assayData,
            pedigree=pedigree,
            ...)
          })

setMethod("geno", "gTrio",
          function(object) {
            assayDataElement(object, "geno")
          })

setMethod("getGeno", "gTrio",
            function(object, ... ){
              if( length(geno(object)) == 0 )
                 stop("No geno data provided.")
              geno.mat.fn(object, ... )
            })

setMethod("pedigree", "gTrio",
            function(object ){
              object@pedigree
            })

setMethod("trios", "gTrio",
            function(object ){
              trios(pedigree(object))
            })

                   #phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   #fatherPhenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   #motherPhenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   #annotation=character(),
                   #featureData=GenomeAnnotatedDataFrameFrom(assayData, annotation, genome=genome),
                   #experimentData=new("MIAME"),
                   #protocolData=phenoData[, integer(0)],
                   #logRRatio=array(NA, dim=c(0, 0, 3)),
                   #BAF=array(NA, dim=c(0,0,3)),
                   #mindist=NULL,
                   #genome=c("hg19", "hg18")
            #.Object@fatherPhenoData <- fatherPhenoData
            #.Object@motherPhenoData <- motherPhenoData
                           #phenoData=phenoData,
                           #fatherPhenoData=fatherPhenoData,
                           #motherPhenoData=motherPhenoData,
                           #featureData=featureData,
                           #experimentData=experimentData,
                           #annotation=annotation,
                           #protocolData=protocolData,
                           #mindist=mindist,
                           #genome=match.arg(genome)
