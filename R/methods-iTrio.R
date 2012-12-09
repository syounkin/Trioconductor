setMethod("initialize", "iTrio",
          function(.Object,
                   assayData = assayDataNew(lrr=lrr, baf = baf),
                   lrr=array(NA, dim=c(0, 0, 3)),
                   baf=array(NA, dim=c(0, 0, 3)),                   
                   pedigree=Pedigree(), ...)
          {
            .Object@pedigree <- pedigree
            callNextMethod(.Object,
            assayData=assayData,
            pedigree=pedigree,
            ...)
          })

