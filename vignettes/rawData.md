    cd ~/trioClasses/vignettes/ && pandoc -o rawData.html rawData.md && cd ~/jhsph/sgy-website/ && ./build-website.sh




```r
library("trioClasses")
library("GWASTools")
library("CleftCNVAssoc")
```






First we create a vector of offspring IDs that we want plotted.

```r
offspring.vec <- as.character(completeTrios(fe.beaty)$id)
```

Now we incorporate it into a GRange object with the ranges in this case being the chr16 region for each of the offspring.

```r
gr <- GRanges(seqnames = rep("chr16", length(offspring.vec)), ranges = IRanges(start = 32404517, 
    end = 32530051), id = offspring.vec)
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
1 0.1662 1.0000 31405309  0.2031  0.0326 1.0000 1.0000    0       0
2 0.3557 1.0000 31405382  0.4870  0.2167 1.0000 1.0000    0       0
3 0.0232 0.9853 31411252  0.0954 -0.0743 0.4786 0.9795    0       1
4 0.1978 0.9981 31428777  0.1213  0.1274 1.0000 1.0000    0       0
5 0.0435 0.0057 31435321 -0.1137  0.1395 0.0071 0.5640    2       2
  geno.ma     x     y  x.fa  y.fa  x.ma  y.ma   snpname
1       0 0.000 1.310 0.000 1.344 0.000 1.194 rs3813007
2       0 0.031 0.734 0.036 0.802 0.025 0.668 rs3813008
3       0 0.194 1.769 1.349 1.270 0.194 1.647 rs4536493
4       0 0.020 0.875 0.005 0.837 0.013 0.836 rs4889545
5       1 1.246 0.058 1.116 0.055 0.763 1.003 rs4477723
```




Plot the logR values for everyone stratified by F,M,O.  Purple is offspring, red is father, and blue is mother.
![plot of chunk logrplot](figures/logrplot1.png) ![plot of chunk logrplot](figures/logrplot2.png) ![plot of chunk logrplot](figures/logrplot3.png) 

Not very informative so we turn to individual trios with an untransmitted deletion. First, we need to find a vector offspring IDs with an untransmitted deletion.  This is a property of the CNVMatrix within the FamilyExperiment object and can be manipulated with the non-exported method TrioAssay.  To begin we first subset the CNVMatrix on the chr16 region.

```r
chr16.gr <- GRanges(seqnames = "chr16", ranges = IRanges(start = 32404517, end = 32530051))
(fe.beaty.chr16 <- fe.beaty[queryHits(findOverlaps(rowData(fe.beaty), chr16.gr))])
```

```
class: FamilyExperiment 
dim: 123 1339 
exptData(0):
assays(1): cnv
rownames(123): comp5899 comp5900 ... comp6020 comp6021
rowData metadata column names(0):
colnames(1339): 11005_01@1008472480 11005_02@1008472482 ...
  18117_02@0070298660 18117_03@0070298657
colData names(1): id
pedigree(2082): famid id fid mid sex dx
complete trios(445):
```

Now with the smaller FE object we can easily construct the trio-states.

```r
trioAssay.chr16 <- trioClasses:::TrioAssay(fe.beaty.chr16, type = "cnv")
trioStates.chr16 <- with(trioAssay.chr16, matrix(paste0(F, M, O), nrow = nrow(O), 
    ncol = ncol(O)))
dimnames(trioStates.chr16) <- dimnames(trioAssay.chr16$O)
head(trioStates.chr16[, 1:5], 10)
```

```
                    comp5899 comp5900 comp5901 comp5902 comp5903
11005_01@1008472480 "000"    "000"    "000"    "000"    "000"   
11021_01@1008472417 "000"    "000"    "000"    "000"    "000"   
11035_01@1008471376 "000"    "000"    "000"    "000"    "000"   
12002_01@1008489061 "000"    "000"    "000"    "000"    "000"   
12004_01@1008489060 "000"    "000"    "000"    "000"    "000"   
12005_01@1008490117 "000"    "000"    "000"    "000"    "000"   
12008_01@1008490140 "010"    "010"    "010"    "010"    "010"   
12014_01@1008490162 "000"    "000"    "000"    "000"    "000"   
12015_01@1008490100 "001"    "001"    "001"    "001"    "001"   
12017_01@1008489083 "000"    "000"    "000"    "000"    "000"   
```

Now we identify trio-cnv pairs with an untransmitted deletion, i.e., trio-states 100, 010, or 110.  (This is not a complete list of trio-states with a non-transmission.)

```r
untrans.mat <- matrix(trioStates.chr16 %in% c("100", "010", "110"), nrow = nrow(trioStates.chr16), 
    ncol = ncol(trioStates.chr16), byrow = FALSE, dimnames = dimnames(trioStates.chr16))
head(untrans.mat[, 1:10], 10)
```

```
                    comp5899 comp5900 comp5901 comp5902 comp5903 comp5904
11005_01@1008472480    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
11021_01@1008472417    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
11035_01@1008471376    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
12002_01@1008489061    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
12004_01@1008489060    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
12005_01@1008490117    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
12008_01@1008490140     TRUE     TRUE     TRUE     TRUE     TRUE     TRUE
12014_01@1008490162    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
12015_01@1008490100    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
12017_01@1008489083    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
                    comp5905 comp5906 comp5907 comp5908
11005_01@1008472480    FALSE    FALSE    FALSE    FALSE
11021_01@1008472417    FALSE    FALSE    FALSE    FALSE
11035_01@1008471376    FALSE    FALSE    FALSE    FALSE
12002_01@1008489061    FALSE    FALSE    FALSE    FALSE
12004_01@1008489060    FALSE    FALSE    FALSE    FALSE
12005_01@1008490117    FALSE    FALSE    FALSE    FALSE
12008_01@1008490140     TRUE     TRUE     TRUE     TRUE
12014_01@1008490162    FALSE    FALSE    FALSE     TRUE
12015_01@1008490100    FALSE    FALSE    FALSE    FALSE
12017_01@1008489083     TRUE     TRUE     TRUE     TRUE
```

And finally we find the IDs of those with more than zero untransmitted deletions.

```r
offspring.chr16 <- rownames(untrans.mat)[which(rowSums(untrans.mat) > 0)]
length(offspring.chr16)
```

```
[1] 140
```

```r
head(offspring.chr16)
```

```
[1] "12008_01@1008490140" "12014_01@1008490162" "12017_01@1008489083"
[4] "12021_01@1008490126" "12024_01@1008490151" "12027_01@1008490157"
```

Now we have what we need to dig deeper into the cause of the under-transmission.
![plot of chunk logrplot2](figures/logrplot21.png) ![plot of chunk logrplot2](figures/logrplot22.png) ![plot of chunk logrplot2](figures/logrplot23.png) ![plot of chunk logrplot2](figures/logrplot24.png) ![plot of chunk logrplot2](figures/logrplot25.png) ![plot of chunk logrplot2](figures/logrplot26.png) ![plot of chunk logrplot2](figures/logrplot27.png) ![plot of chunk logrplot2](figures/logrplot28.png) ![plot of chunk logrplot2](figures/logrplot29.png) ![plot of chunk logrplot2](figures/logrplot210.png) ![plot of chunk logrplot2](figures/logrplot211.png) ![plot of chunk logrplot2](figures/logrplot212.png) ![plot of chunk logrplot2](figures/logrplot213.png) ![plot of chunk logrplot2](figures/logrplot214.png) ![plot of chunk logrplot2](figures/logrplot215.png) ![plot of chunk logrplot2](figures/logrplot216.png) ![plot of chunk logrplot2](figures/logrplot217.png) ![plot of chunk logrplot2](figures/logrplot218.png) ![plot of chunk logrplot2](figures/logrplot219.png) ![plot of chunk logrplot2](figures/logrplot220.png) ![plot of chunk logrplot2](figures/logrplot221.png) ![plot of chunk logrplot2](figures/logrplot222.png) ![plot of chunk logrplot2](figures/logrplot223.png) ![plot of chunk logrplot2](figures/logrplot224.png) ![plot of chunk logrplot2](figures/logrplot225.png) ![plot of chunk logrplot2](figures/logrplot226.png) ![plot of chunk logrplot2](figures/logrplot227.png) ![plot of chunk logrplot2](figures/logrplot228.png) ![plot of chunk logrplot2](figures/logrplot229.png) ![plot of chunk logrplot2](figures/logrplot230.png) ![plot of chunk logrplot2](figures/logrplot231.png) ![plot of chunk logrplot2](figures/logrplot232.png) ![plot of chunk logrplot2](figures/logrplot233.png) ![plot of chunk logrplot2](figures/logrplot234.png) ![plot of chunk logrplot2](figures/logrplot235.png) ![plot of chunk logrplot2](figures/logrplot236.png) ![plot of chunk logrplot2](figures/logrplot237.png) ![plot of chunk logrplot2](figures/logrplot238.png) ![plot of chunk logrplot2](figures/logrplot239.png) ![plot of chunk logrplot2](figures/logrplot240.png) ![plot of chunk logrplot2](figures/logrplot241.png) ![plot of chunk logrplot2](figures/logrplot242.png) ![plot of chunk logrplot2](figures/logrplot243.png) ![plot of chunk logrplot2](figures/logrplot244.png) ![plot of chunk logrplot2](figures/logrplot245.png) ![plot of chunk logrplot2](figures/logrplot246.png) ![plot of chunk logrplot2](figures/logrplot247.png) ![plot of chunk logrplot2](figures/logrplot248.png) ![plot of chunk logrplot2](figures/logrplot249.png) ![plot of chunk logrplot2](figures/logrplot250.png) ![plot of chunk logrplot2](figures/logrplot251.png) ![plot of chunk logrplot2](figures/logrplot252.png) ![plot of chunk logrplot2](figures/logrplot253.png) ![plot of chunk logrplot2](figures/logrplot254.png) ![plot of chunk logrplot2](figures/logrplot255.png) ![plot of chunk logrplot2](figures/logrplot256.png) ![plot of chunk logrplot2](figures/logrplot257.png) ![plot of chunk logrplot2](figures/logrplot258.png) ![plot of chunk logrplot2](figures/logrplot259.png) ![plot of chunk logrplot2](figures/logrplot260.png) ![plot of chunk logrplot2](figures/logrplot261.png) ![plot of chunk logrplot2](figures/logrplot262.png) ![plot of chunk logrplot2](figures/logrplot263.png) ![plot of chunk logrplot2](figures/logrplot264.png) ![plot of chunk logrplot2](figures/logrplot265.png) ![plot of chunk logrplot2](figures/logrplot266.png) ![plot of chunk logrplot2](figures/logrplot267.png) ![plot of chunk logrplot2](figures/logrplot268.png) ![plot of chunk logrplot2](figures/logrplot269.png) ![plot of chunk logrplot2](figures/logrplot270.png) ![plot of chunk logrplot2](figures/logrplot271.png) ![plot of chunk logrplot2](figures/logrplot272.png) ![plot of chunk logrplot2](figures/logrplot273.png) ![plot of chunk logrplot2](figures/logrplot274.png) ![plot of chunk logrplot2](figures/logrplot275.png) ![plot of chunk logrplot2](figures/logrplot276.png) ![plot of chunk logrplot2](figures/logrplot277.png) ![plot of chunk logrplot2](figures/logrplot278.png) ![plot of chunk logrplot2](figures/logrplot279.png) ![plot of chunk logrplot2](figures/logrplot280.png) ![plot of chunk logrplot2](figures/logrplot281.png) ![plot of chunk logrplot2](figures/logrplot282.png) ![plot of chunk logrplot2](figures/logrplot283.png) ![plot of chunk logrplot2](figures/logrplot284.png) ![plot of chunk logrplot2](figures/logrplot285.png) ![plot of chunk logrplot2](figures/logrplot286.png) ![plot of chunk logrplot2](figures/logrplot287.png) ![plot of chunk logrplot2](figures/logrplot288.png) ![plot of chunk logrplot2](figures/logrplot289.png) ![plot of chunk logrplot2](figures/logrplot290.png) ![plot of chunk logrplot2](figures/logrplot291.png) ![plot of chunk logrplot2](figures/logrplot292.png) ![plot of chunk logrplot2](figures/logrplot293.png) ![plot of chunk logrplot2](figures/logrplot294.png) ![plot of chunk logrplot2](figures/logrplot295.png) ![plot of chunk logrplot2](figures/logrplot296.png) ![plot of chunk logrplot2](figures/logrplot297.png) ![plot of chunk logrplot2](figures/logrplot298.png) ![plot of chunk logrplot2](figures/logrplot299.png) ![plot of chunk logrplot2](figures/logrplot2100.png) ![plot of chunk logrplot2](figures/logrplot2101.png) ![plot of chunk logrplot2](figures/logrplot2102.png) ![plot of chunk logrplot2](figures/logrplot2103.png) ![plot of chunk logrplot2](figures/logrplot2104.png) ![plot of chunk logrplot2](figures/logrplot2105.png) ![plot of chunk logrplot2](figures/logrplot2106.png) ![plot of chunk logrplot2](figures/logrplot2107.png) ![plot of chunk logrplot2](figures/logrplot2108.png) ![plot of chunk logrplot2](figures/logrplot2109.png) ![plot of chunk logrplot2](figures/logrplot2110.png) ![plot of chunk logrplot2](figures/logrplot2111.png) ![plot of chunk logrplot2](figures/logrplot2112.png) ![plot of chunk logrplot2](figures/logrplot2113.png) ![plot of chunk logrplot2](figures/logrplot2114.png) ![plot of chunk logrplot2](figures/logrplot2115.png) ![plot of chunk logrplot2](figures/logrplot2116.png) ![plot of chunk logrplot2](figures/logrplot2117.png) ![plot of chunk logrplot2](figures/logrplot2118.png) ![plot of chunk logrplot2](figures/logrplot2119.png) ![plot of chunk logrplot2](figures/logrplot2120.png) ![plot of chunk logrplot2](figures/logrplot2121.png) ![plot of chunk logrplot2](figures/logrplot2122.png) ![plot of chunk logrplot2](figures/logrplot2123.png) ![plot of chunk logrplot2](figures/logrplot2124.png) ![plot of chunk logrplot2](figures/logrplot2125.png) ![plot of chunk logrplot2](figures/logrplot2126.png) ![plot of chunk logrplot2](figures/logrplot2127.png) ![plot of chunk logrplot2](figures/logrplot2128.png) ![plot of chunk logrplot2](figures/logrplot2129.png) ![plot of chunk logrplot2](figures/logrplot2130.png) ![plot of chunk logrplot2](figures/logrplot2131.png) ![plot of chunk logrplot2](figures/logrplot2132.png) ![plot of chunk logrplot2](figures/logrplot2133.png) ![plot of chunk logrplot2](figures/logrplot2134.png) ![plot of chunk logrplot2](figures/logrplot2135.png) ![plot of chunk logrplot2](figures/logrplot2136.png) ![plot of chunk logrplot2](figures/logrplot2137.png) ![plot of chunk logrplot2](figures/logrplot2138.png) ![plot of chunk logrplot2](figures/logrplot2139.png) ![plot of chunk logrplot2](figures/logrplot2140.png) 

Now for the transmitted deletions in the chr16 region.

```r
trans.mat <- matrix(trioStates.chr16 %in% c("101", "011", "112"), nrow = nrow(trioStates.chr16), 
    ncol = ncol(trioStates.chr16), byrow = FALSE, dimnames = dimnames(trioStates.chr16))
```


```r
offspring.chr16 <- rownames(trans.mat)[which(rowSums(trans.mat) > 0)]
length(offspring.chr16)
```

```
[1] 69
```

```r
head(offspring.chr16)
```

```
[1] "12008_01@1008490140" "12014_01@1008490162" "12054_01@1008494951"
[4] "12062_01@0067868215" "12064_01@0067868240" "12071_01@0067868170"
```

![plot of chunk logrplot3](figures/logrplot31.png) ![plot of chunk logrplot3](figures/logrplot32.png) ![plot of chunk logrplot3](figures/logrplot33.png) ![plot of chunk logrplot3](figures/logrplot34.png) ![plot of chunk logrplot3](figures/logrplot35.png) ![plot of chunk logrplot3](figures/logrplot36.png) ![plot of chunk logrplot3](figures/logrplot37.png) ![plot of chunk logrplot3](figures/logrplot38.png) ![plot of chunk logrplot3](figures/logrplot39.png) ![plot of chunk logrplot3](figures/logrplot310.png) ![plot of chunk logrplot3](figures/logrplot311.png) ![plot of chunk logrplot3](figures/logrplot312.png) ![plot of chunk logrplot3](figures/logrplot313.png) ![plot of chunk logrplot3](figures/logrplot314.png) ![plot of chunk logrplot3](figures/logrplot315.png) ![plot of chunk logrplot3](figures/logrplot316.png) ![plot of chunk logrplot3](figures/logrplot317.png) ![plot of chunk logrplot3](figures/logrplot318.png) ![plot of chunk logrplot3](figures/logrplot319.png) ![plot of chunk logrplot3](figures/logrplot320.png) ![plot of chunk logrplot3](figures/logrplot321.png) ![plot of chunk logrplot3](figures/logrplot322.png) ![plot of chunk logrplot3](figures/logrplot323.png) ![plot of chunk logrplot3](figures/logrplot324.png) ![plot of chunk logrplot3](figures/logrplot325.png) ![plot of chunk logrplot3](figures/logrplot326.png) ![plot of chunk logrplot3](figures/logrplot327.png) ![plot of chunk logrplot3](figures/logrplot328.png) ![plot of chunk logrplot3](figures/logrplot329.png) ![plot of chunk logrplot3](figures/logrplot330.png) ![plot of chunk logrplot3](figures/logrplot331.png) ![plot of chunk logrplot3](figures/logrplot332.png) ![plot of chunk logrplot3](figures/logrplot333.png) ![plot of chunk logrplot3](figures/logrplot334.png) ![plot of chunk logrplot3](figures/logrplot335.png) ![plot of chunk logrplot3](figures/logrplot336.png) ![plot of chunk logrplot3](figures/logrplot337.png) ![plot of chunk logrplot3](figures/logrplot338.png) ![plot of chunk logrplot3](figures/logrplot339.png) ![plot of chunk logrplot3](figures/logrplot340.png) ![plot of chunk logrplot3](figures/logrplot341.png) ![plot of chunk logrplot3](figures/logrplot342.png) ![plot of chunk logrplot3](figures/logrplot343.png) ![plot of chunk logrplot3](figures/logrplot344.png) ![plot of chunk logrplot3](figures/logrplot345.png) ![plot of chunk logrplot3](figures/logrplot346.png) ![plot of chunk logrplot3](figures/logrplot347.png) ![plot of chunk logrplot3](figures/logrplot348.png) ![plot of chunk logrplot3](figures/logrplot349.png) ![plot of chunk logrplot3](figures/logrplot350.png) ![plot of chunk logrplot3](figures/logrplot351.png) ![plot of chunk logrplot3](figures/logrplot352.png) ![plot of chunk logrplot3](figures/logrplot353.png) ![plot of chunk logrplot3](figures/logrplot354.png) ![plot of chunk logrplot3](figures/logrplot355.png) ![plot of chunk logrplot3](figures/logrplot356.png) ![plot of chunk logrplot3](figures/logrplot357.png) ![plot of chunk logrplot3](figures/logrplot358.png) ![plot of chunk logrplot3](figures/logrplot359.png) ![plot of chunk logrplot3](figures/logrplot360.png) ![plot of chunk logrplot3](figures/logrplot361.png) ![plot of chunk logrplot3](figures/logrplot362.png) ![plot of chunk logrplot3](figures/logrplot363.png) ![plot of chunk logrplot3](figures/logrplot364.png) ![plot of chunk logrplot3](figures/logrplot365.png) ![plot of chunk logrplot3](figures/logrplot366.png) ![plot of chunk logrplot3](figures/logrplot367.png) ![plot of chunk logrplot3](figures/logrplot368.png) ![plot of chunk logrplot3](figures/logrplot369.png) 

