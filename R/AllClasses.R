setClass("PedClass", contains = "DataFrame", prototype = DataFrame( famid=factor(), id=factor(), fid=factor(), mid=factor(), sex=factor(), dx=factor()))
setClass("FamilyExperiment", contains="SummarizedExperiment", representation(pedigree="PedClass") )
