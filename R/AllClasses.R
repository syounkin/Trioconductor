setClass("PedClass", contains = "DataFrame", prototype = DataFrame( famid=paste0("fam",rep(1:3,3)), id=paste0("pseudo",1:9), fid=paste0("pseudo",1:9), mid=paste0("pseudo",1:9), sex=paste0("pseudo",1:9), dx=paste0("pseudo",1:9)))

setClass("SNPTrioExperiment", contains="SummarizedExperiment", representation(pedigree="PedClass") )
