


```r
library("trioClasses")
library("GWASTools")
library("CleftCNVAssoc")
```






First we create a vector of offspring IDs that we want plotted.

```r
offspring.vec <- as.character(completeTrios(fe.beaty)$id)
```

Now we incorporate it into a GRange object with the ranges in this case being the chr15 region for each of the offspring.

```r
gr <- GRanges(seqnames = rep("chr15", length(offspring.vec)), ranges = IRanges(start = 19768826, 
    end = 19982036), id = offspring.vec)
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
1 -0.2163 0.0000 18788683 -0.0209  0.2160 0.0000 0.0616   NA      NA
2 -0.1305 0.4200 18822301  0.1224  0.1599 0.4511 0.4183   NA      NA
3 -0.2070 0.9655 18832819  0.2894  0.1535 0.5595 0.5668   NA      NA
4 -0.3455 0.0045 18850028 -0.0927 -0.0213 0.0107 0.0177   NA      NA
5 -0.0958 0.9988 18851201  0.0892 -0.0873 0.9916 1.0000    0       0
  geno.ma     x     y  x.fa  y.fa  x.ma  y.ma    snpname
1      NA 1.063 0.003 1.218 0.003 1.290 0.149  rs4114751
2      NA 0.611 1.198 0.636 1.520 0.752 1.461  rs4414456
3      NA 0.048 0.578 0.194 0.689 0.175 0.628  rs1813365
4      NA 1.300 0.061 1.537 0.085 1.600 0.104  rs8041088
5       0 0.046 1.128 0.064 1.271 0.042 1.140 rs28690063
```




Plot the logR values for everyone stratified by F,M,O.  Purple is offspring, red is father, and blue is mother.


Not very informative so we turn to individual trios with an untransmitted deletion. First, we need to find a vector offspring IDs with an untransmitted deletion.  This is a property of the CNVMatrix within the FamilyExperiment object and can be manipulated with the non-exported method TrioAssay.  To begin we first subset the CNVMatrix on the chr15 region.

```r
chr15.gr <- GRanges(seqnames = "chr15", ranges = IRanges(start = 19768826, end = 19982036))

(fe.beaty.chr15 <- fe.beaty[queryHits(findOverlaps(rowData(fe.beaty), chr15.gr))])
```

```
class: FamilyExperiment 
dim: 21 1339 
exptData(0):
assays(1): cnv
rownames(21): comp5379 comp5380 ... comp5398 comp5399
rowData metadata column names(0):
colnames(1339): 11005_01@1008472480 11005_02@1008472482 ...
  18117_02@0070298660 18117_03@0070298657
colData names(1): id
pedigree(2082): famid id fid mid sex dx
complete trios(445):
```

Now with the smaller FE object we can easily construct the trio-states.

```r
trioAssay.chr15 <- trioClasses:::TrioAssay(fe.beaty.chr15, type = "cnv")
trioStates.chr15 <- with(trioAssay.chr15, matrix(paste0(F, M, O), nrow = nrow(O), 
    ncol = ncol(O)))
dimnames(trioStates.chr15) <- dimnames(trioAssay.chr15$O)
head(trioStates.chr15[, 1:5], 10)
```

```
                    comp5379 comp5380 comp5381 comp5382 comp5383
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
untrans.mat <- matrix(trioStates.chr15 %in% c("100", "010", "110"), nrow = nrow(trioStates.chr15), 
    ncol = ncol(trioStates.chr15), byrow = FALSE, dimnames = dimnames(trioStates.chr15))
head(untrans.mat[, 1:10], 10)
```

```
                    comp5379 comp5380 comp5381 comp5382 comp5383 comp5384
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
                    comp5385 comp5386 comp5387 comp5388
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
offspring.chr15 <- rownames(untrans.mat)[which(rowSums(untrans.mat) > 0)]
length(offspring.chr15)
```

```
[1] 31
```

```r
head(offspring.chr15)
```

```
[1] "12054_01@1008494951" "12081_01@0067868174" "12093_01@0067866270"
[4] "12102_01@0067866338" "12158_01@0067866940" "12173_01@0067866951"
```

![plot of chunk logrplot2-chr15](figures/logrplot2-chr151.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr152.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr153.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr154.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr155.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr156.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr157.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr158.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr159.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1510.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1511.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1512.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1513.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1514.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1515.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1516.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1517.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1518.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1519.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1520.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1521.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1522.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1523.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1524.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1525.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1526.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1527.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1528.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1529.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1530.png) ![plot of chunk logrplot2-chr15](figures/logrplot2-chr1531.png) 

Now for the transmitted deletions in the chr15 region.

```r
trans.mat <- matrix(trioStates.chr15 %in% c("101", "011", "112"), nrow = nrow(trioStates.chr15), 
    ncol = ncol(trioStates.chr15), byrow = FALSE, dimnames = dimnames(trioStates.chr15))
```


```r
offspring.chr15 <- rownames(trans.mat)[which(rowSums(trans.mat) > 0)]
length(offspring.chr15)
```

```
[1] 31
```

```r
head(offspring.chr15)
```

```
[1] "12059_01@1008494932" "12066_01@0067868260" "12112_01@0067866334"
[4] "12137_01@0067867189" "12151_01@0067866460" "12157_01@0067866969"
```

![plot of chunk logrplot3-chr15](figures/logrplot3-chr151.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr152.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr153.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr154.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr155.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr156.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr157.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr158.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr159.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1510.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1511.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1512.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1513.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1514.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1515.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1516.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1517.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1518.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1519.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1520.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1521.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1522.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1523.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1524.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1525.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1526.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1527.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1528.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1529.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1530.png) ![plot of chunk logrplot3-chr15](figures/logrplot3-chr1531.png) 

