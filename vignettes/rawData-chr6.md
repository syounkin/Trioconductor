


```r
library("trioClasses")
library("GWASTools")
library("CleftCNVAssoc")
```






First we create a vector of offspring IDs that we want plotted.

```r
offspring.vec <- as.character(completeTrios(fe.beaty)$id)
```

Now we incorporate it into a GRange object with the ranges in this case being the chr6 region for each of the offspring.

```r
gr <- GRanges(seqnames = rep("chr6", length(offspring.vec)), ranges = IRanges(start = 32611466, 
    end = 32643872), id = offspring.vec)
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
     logr    baf      pos logr.fa logr.ma baf.fa baf.ma geno geno.fa
1 -0.0035 0.0039 31611810  0.0561 -0.1128 0.0061 0.0059    2       2
2 -0.1510 1.0000 31612525 -0.0740 -0.1015 1.0000 0.9769    0       0
3 -0.1297 0.0189 31612566  0.2524 -0.0882 0.0030 0.0070    2       2
4  0.0263 0.0000 31613459  0.1762  0.1724 0.0009 0.0000    2       2
5  0.0552 0.9987 31613866  0.0433  0.0543 0.9819 0.9955    0       0
  geno.ma     x     y  x.fa  y.fa  x.ma  y.ma    snpname
1       2 1.994 0.024 2.075 0.031 1.846 0.027 rs17200712
2       0 0.059 1.179 0.076 1.233 0.113 1.179 rs10456396
3       2 0.513 0.041 0.682 0.042 0.536 0.035  rs3130057
4       2 0.696 0.035 0.770 0.042 0.775 0.033  rs2734583
5       0 0.030 1.414 0.060 1.384 0.036 1.409  rs3130058
```




Plot the logR values for everyone stratified by F,M,O.  Purple is offspring, red is father, and blue is mother.


Not very informative so we turn to individual trios with an untransmitted deletion. First, we need to find a vector offspring IDs with an untransmitted deletion.  This is a property of the CNVMatrix within the FamilyExperiment object and can be manipulated with the non-exported method TrioAssay.  To begin we first subset the CNVMatrix on the chr6 region.

```r
chr6.gr <- GRanges(seqnames = "chr6", ranges = IRanges(start = 32611466, end = 32643872))

(fe.beaty.chr6 <- fe.beaty[queryHits(findOverlaps(rowData(fe.beaty), chr6.gr))])
```

```
class: FamilyExperiment 
dim: 43 1339 
exptData(0):
assays(1): cnv
rownames(43): comp2326 comp2327 ... comp2367 comp2368
rowData metadata column names(0):
colnames(1339): 11005_01@1008472480 11005_02@1008472482 ...
  18117_02@0070298660 18117_03@0070298657
colData names(1): id
pedigree(2082): famid id fid mid sex dx
complete trios(445):
```

Now with the smaller FE object we can easily construct the trio-states.

```r
trioAssay.chr6 <- trioClasses:::TrioAssay(fe.beaty.chr6, type = "cnv")
trioStates.chr6 <- with(trioAssay.chr6, matrix(paste0(F, M, O), nrow = nrow(O), 
    ncol = ncol(O)))
dimnames(trioStates.chr6) <- dimnames(trioAssay.chr6$O)
head(trioStates.chr6[, 1:5], 10)
```

```
                    comp2326 comp2327 comp2328 comp2329 comp2330
11005_01@1008472480 "000"    "000"    "000"    "000"    "000"   
11021_01@1008472417 "000"    "000"    "000"    "000"    "000"   
11035_01@1008471376 "000"    "000"    "000"    "000"    "000"   
12002_01@1008489061 "000"    "000"    "000"    "000"    "000"   
12004_01@1008489060 "000"    "000"    "000"    "000"    "000"   
12005_01@1008490117 "000"    "000"    "000"    "000"    "000"   
12008_01@1008490140 "000"    "000"    "000"    "000"    "000"   
12014_01@1008490162 "000"    "000"    "000"    "000"    "000"   
12015_01@1008490100 "000"    "000"    "000"    "000"    "000"   
12017_01@1008489083 "000"    "000"    "000"    "000"    "000"   
```

Now we identify trio-cnv pairs with an untransmitted deletion, i.e., trio-states 100, 010, or 110.  (This is not a complete list of trio-states with a non-transmission.)

```r
untrans.mat <- matrix(trioStates.chr6 %in% c("100", "010", "110"), nrow = nrow(trioStates.chr6), 
    ncol = ncol(trioStates.chr6), byrow = FALSE, dimnames = dimnames(trioStates.chr6))
head(untrans.mat[, 1:10], 10)
```

```
                    comp2326 comp2327 comp2328 comp2329 comp2330 comp2331
11005_01@1008472480    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
11021_01@1008472417    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
11035_01@1008471376    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
12002_01@1008489061    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
12004_01@1008489060    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
12005_01@1008490117    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
12008_01@1008490140    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
12014_01@1008490162    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
12015_01@1008490100    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
12017_01@1008489083    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
                    comp2332 comp2333 comp2334 comp2335
11005_01@1008472480    FALSE    FALSE    FALSE    FALSE
11021_01@1008472417    FALSE    FALSE    FALSE    FALSE
11035_01@1008471376    FALSE    FALSE    FALSE    FALSE
12002_01@1008489061    FALSE    FALSE    FALSE    FALSE
12004_01@1008489060    FALSE    FALSE    FALSE    FALSE
12005_01@1008490117    FALSE    FALSE    FALSE    FALSE
12008_01@1008490140    FALSE    FALSE    FALSE    FALSE
12014_01@1008490162    FALSE    FALSE    FALSE    FALSE
12015_01@1008490100    FALSE    FALSE    FALSE    FALSE
12017_01@1008489083    FALSE    FALSE    FALSE    FALSE
```

And finally we find the IDs of those with more than zero untransmitted deletions.

```r
offspring.chr6 <- rownames(untrans.mat)[which(rowSums(untrans.mat) > 0)]
length(offspring.chr6)
```

```
[1] 65
```

```r
head(offspring.chr6)
```

```
[1] "12015_01@1008490100" "12017_01@1008489083" "12039_01@1008494891"
[4] "12042_01@1008494899" "12047_01@1008494888" "12059_01@1008494932"
```

Now we have what we need to dig deeper into the cause of the under-transmission.
![plot of chunk logrplot2-chr6](figures/logrplot2-chr61.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr62.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr63.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr64.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr65.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr66.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr67.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr68.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr69.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr610.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr611.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr612.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr613.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr614.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr615.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr616.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr617.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr618.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr619.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr620.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr621.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr622.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr623.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr624.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr625.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr626.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr627.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr628.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr629.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr630.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr631.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr632.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr633.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr634.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr635.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr636.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr637.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr638.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr639.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr640.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr641.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr642.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr643.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr644.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr645.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr646.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr647.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr648.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr649.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr650.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr651.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr652.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr653.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr654.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr655.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr656.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr657.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr658.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr659.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr660.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr661.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr662.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr663.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr664.png) ![plot of chunk logrplot2-chr6](figures/logrplot2-chr665.png) 

Now for the transmitted deletions in the chr6 region.

```r
trans.mat <- matrix(trioStates.chr6 %in% c("101", "011", "112"), nrow = nrow(trioStates.chr6), 
    ncol = ncol(trioStates.chr6), byrow = FALSE, dimnames = dimnames(trioStates.chr6))
```


```r
offspring.chr6 <- rownames(trans.mat)[which(rowSums(trans.mat) > 0)]
length(offspring.chr6)
```

```
[1] 45
```

```r
head(offspring.chr6)
```

```
[1] "12039_01@1008494891" "12052_01@1008494879" "12054_01@1008494951"
[4] "12080_01@0067868222" "12089_01@0067868251" "12095_01@0067866272"
```

![plot of chunk logrplot3-chr6](figures/logrplot3-chr61.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr62.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr63.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr64.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr65.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr66.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr67.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr68.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr69.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr610.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr611.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr612.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr613.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr614.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr615.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr616.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr617.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr618.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr619.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr620.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr621.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr622.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr623.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr624.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr625.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr626.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr627.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr628.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr629.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr630.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr631.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr632.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr633.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr634.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr635.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr636.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr637.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr638.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr639.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr640.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr641.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr642.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr643.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr644.png) ![plot of chunk logrplot3-chr6](figures/logrplot3-chr645.png) 

