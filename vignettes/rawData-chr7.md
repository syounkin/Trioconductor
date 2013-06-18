


```r
library("trioClasses")
library("GWASTools")
library("CleftCNVAssoc")
```






First we create a vector of offspring IDs that we want plotted.

```r
offspring.vec <- as.character(completeTrios(fe.beaty)$id)
```

Now we incorporate it into a GRange object with the ranges in this case being the chr7 region for each of the offspring.

```r
gr <- GRanges(seqnames = rep("chr7", length(offspring.vec)), ranges = IRanges(start = 141419097, 
    end = 141441259), id = offspring.vec)
```

Now we apply a function from CleftCNVAssoc to retrieve data from GWASTools objects.  These data exist only on an encrypted hard drive and enigma.

```r
raw.df.list <- getRaw(gr + 1e+06, intensfile = intensfile, snpAnnot = beaty_snpAnnot, 
    scan.id = scan.ids, fa.id = fa.id, ma.id = ma.id, genofile = genofile, xyfile = xyfile)
```

Here are the data for the first offspring's first 5 markers.

```r
head(raw.df.list[[1]], 5)
```

```
     logr    baf       pos logr.fa logr.ma baf.fa baf.ma geno geno.fa
1 -0.1648 0.5189 140422294 -0.2148 -0.0340 0.0011 1.0000    1       2
2 -0.0777 0.0069 140423362 -0.0257 -0.0588 0.0041 0.0040   NA      NA
3  0.0018 0.9972 140426209  0.1944  0.0264 0.5048 1.0000    0       1
4  0.0450 0.5078 140426359  0.0468  0.1053 0.0000 0.4995    1       2
5  0.0962 0.9957 140433319 -0.0501 -0.0142 1.0000 1.0000   NA      NA
  geno.ma     x     y  x.fa  y.fa  x.ma  y.ma     snpname
1       0 0.762 0.839 1.261 0.016 0.054 1.377    rs498933
2      NA 0.646 0.032 0.673 0.031 0.657 0.030 cnvi0014408
3       0 0.012 0.907 0.585 0.641 0.000 0.932    rs557962
4       1 0.982 0.869 1.600 0.012 1.041 0.896   rs4726259
5      NA 0.018 0.611 0.008 0.561 0.010 0.573 cnvi0014398
```




Plot the logR values for everyone stratified by F,M,O.  Purple is offspring, red is father, and blue is mother.


Not very informative so we turn to individual trios with an untransmitted deletion. First, we need to find a vector offspring IDs with an untransmitted deletion.  This is a property of the CNVMatrix within the FamilyExperiment object and can be manipulated with the non-exported method TrioAssay.  To begin we first subset the CNVMatrix on the chr7 region.

```r
chr7.gr <- GRanges(seqnames = "chr7", ranges = IRanges(start = 141419097, end = 141441259))

(fe.beaty.chr7 <- fe.beaty[queryHits(findOverlaps(rowData(fe.beaty), chr7.gr))])
```

```
class: FamilyExperiment 
dim: 8 1339 
exptData(0):
assays(1): cnv
rownames(8): comp3086 comp3087 ... comp3092 comp3093
rowData metadata column names(0):
colnames(1339): 11005_01@1008472480 11005_02@1008472482 ...
  18117_02@0070298660 18117_03@0070298657
colData names(1): id
pedigree(2082): famid id fid mid sex dx
complete trios(445):
```

Now with the smaller FE object we can easily construct the trio-states.

```r
trioAssay.chr7 <- trioClasses:::TrioAssay(fe.beaty.chr7, type = "cnv")
trioStates.chr7 <- with(trioAssay.chr7, matrix(paste0(F, M, O), nrow = nrow(O), 
    ncol = ncol(O)))
dimnames(trioStates.chr7) <- dimnames(trioAssay.chr7$O)
# head(trioStates.chr7[,1:5],10)
```

Now we identify trio-cnv pairs with an untransmitted deletion, i.e., trio-states 100, 010, or 110.  (This is not a complete list of trio-states with a non-transmission.)

```r
untrans.mat <- matrix(trioStates.chr7 %in% c("100", "010", "110"), nrow = nrow(trioStates.chr7), 
    ncol = ncol(trioStates.chr7), byrow = FALSE, dimnames = dimnames(trioStates.chr7))
head(untrans.mat[, 1:10], 10)
```

```
Error: error in evaluating the argument 'x' in selecting a method for
function 'head': Error in untrans.mat[, 1:10] (from <text>#3) : subscript
out of bounds
```

And finally we find the IDs of those with more than zero untransmitted deletions.

```r
offspring.chr7 <- rownames(untrans.mat)[which(rowSums(untrans.mat) > 0)]
length(offspring.chr7)
```

```
[1] 44
```

```r
head(offspring.chr7)
```

```
[1] "12004_01@1008489060" "12026_01@1008490128" "12036_01@1008494873"
[4] "12043_01@1008494947" "12071_01@0067868170" "12074_01@0067868244"
```

Now we have what we need to dig deeper into the cause of the under-transmission.
![plot of chunk logrplot2-chr7](figures/logrplot2-chr71.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr72.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr73.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr74.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr75.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr76.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr77.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr78.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr79.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr710.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr711.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr712.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr713.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr714.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr715.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr716.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr717.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr718.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr719.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr720.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr721.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr722.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr723.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr724.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr725.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr726.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr727.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr728.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr729.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr730.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr731.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr732.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr733.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr734.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr735.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr736.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr737.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr738.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr739.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr740.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr741.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr742.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr743.png) ![plot of chunk logrplot2-chr7](figures/logrplot2-chr744.png) 

Now for the transmitted deletions in the chr7 region.

```r
trans.mat <- matrix(trioStates.chr7 %in% c("101", "011", "112"), nrow = nrow(trioStates.chr7), 
    ncol = ncol(trioStates.chr7), byrow = FALSE, dimnames = dimnames(trioStates.chr7))
```


```r
offspring.chr7 <- rownames(trans.mat)[which(rowSums(trans.mat) > 0)]
length(offspring.chr7)
```

```
[1] 62
```

```r
head(offspring.chr7)
```

```
[1] "12037_01@1008494916" "12040_01@1008494922" "12059_01@1008494932"
[4] "12091_01@0067866342" "12099_01@0067866249" "12107_01@0067866336"
```

![plot of chunk logrplot3-chr7](figures/logrplot3-chr71.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr72.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr73.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr74.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr75.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr76.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr77.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr78.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr79.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr710.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr711.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr712.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr713.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr714.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr715.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr716.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr717.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr718.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr719.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr720.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr721.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr722.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr723.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr724.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr725.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr726.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr727.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr728.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr729.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr730.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr731.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr732.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr733.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr734.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr735.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr736.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr737.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr738.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr739.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr740.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr741.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr742.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr743.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr744.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr745.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr746.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr747.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr748.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr749.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr750.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr751.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr752.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr753.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr754.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr755.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr756.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr757.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr758.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr759.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr760.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr761.png) ![plot of chunk logrplot3-chr7](figures/logrplot3-chr762.png) 

