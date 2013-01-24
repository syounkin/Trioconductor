setMethod("initialize", "SNPExperiment",
          function(.Object,
                   assays = NA,
                   colData = NA,
                   rowData = NA,
                   Pedigree = NA,
                   ...){
            .Object@Pedigree <- Pedigree
            callNextMethod(.Object, assays=assays, colData = colData, rowData = rowData, ...)
          })
