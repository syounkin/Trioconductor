

```r
opts_chunk$set(fig.width = 5, fig.height = 5, width = 200, continue = " ")
library("trioClasses")
library("Gviz")
library("TxDb.Hsapiens.UCSC.hg18.knownGene")
data("fe", package = "trioClasses")
```

Create trio-states

```r
trioAssay.beaty <- trioClasses:::TrioAssay(fe.beaty, type = "cnv")
trioStates.beaty <- with(trioAssay.beaty, matrix(paste0(F, M, O), nrow = nrow(O), 
    ncol = ncol(O)))
dimnames(trioStates.beaty) <- dimnames(trioAssay.beaty$O)
trioAssay.pitt <- trioClasses:::TrioAssay(fe.pitt, type = "cnv")
trioStates.pitt <- with(trioAssay.pitt, matrix(paste0(F, M, O), nrow = nrow(O), 
    ncol = ncol(O)))
dimnames(trioStates.pitt) <- dimnames(trioAssay.pitt$O)
```

Tabulate trio-states component-wise.

```r
table.list.beaty <- apply(trioStates.beaty, 2, "table")
table.list.pitt <- apply(trioStates.pitt, 2, "table")
```

Count transmission events and  compute p-values component-wise.

```r
TU.mat.beaty <- matrix(unlist(lapply(table.list.beaty, trioClasses:::CountTU)), 
    nrow = length(table.list.beaty), ncol = 2, byrow = TRUE)
TU.mat.pitt <- matrix(unlist(lapply(table.list.pitt, trioClasses:::CountTU)), 
    nrow = length(table.list.pitt), ncol = 2, byrow = TRUE)
TU.mat <- cbind(TU.mat.beaty, TU.mat.pitt)
testable <- which((rowSums(TU.mat[, 1:2]) >= 25) & (rowSums(TU.mat[, 3:4]) >= 
    25))
TU.mat <- TU.mat[testable, ]
rownames(TU.mat) <- names(table.list.beaty)[testable]
colnames(TU.mat) <- c("T.case", "U.case", "T.con", "U.con")
DF <- DataFrame(rowData(fe.beaty)[testable], TU.mat)
colnames(DF) <- c("grange", colnames(TU.mat))
```


```r
hist(trans.vec <- rowSums(TU.mat[, c(1, 3)])/rowSums(TU.mat), breaks = 20)
```

![plot of chunk hist](figure/hist.png) 


```r
fish.list <- apply(TU.mat, 1, trioClasses:::TU.fish)
p.vec <- unlist(lapply(fish.list, function(obj) return(obj$p.value)))
DF <- DataFrame(DF, p.vec, trans.vec)
```

Now we look at regions.

```r
regions.gr <- reduce(DF$grange)
p.min.DF <- trioClasses:::f.cmp(DF, "p.vec", min, na.rm = TRUE)
n.DF <- trioClasses:::f.cmp(DF, "p.vec", function(vec) {
    sum(!is.na(vec), na.rm = TRUE)
})
p.median.DF <- trioClasses:::f.cmp(DF, "p.vec", median, na.rm = TRUE)
trans.median.DF <- trioClasses:::f.cmp(DF, "trans.vec", median, na.rm = TRUE)
meta <- values(regions.gr)
meta <- DataFrame(meta, p.min = p.min.DF$value, p.median = p.median.DF$value, 
    trans.median = trans.median.DF$value, n.cmp = n.DF$value)
values(regions.gr) <- meta
```


```r
head(as(regions.gr[order(values(regions.gr)$p.min)], "data.frame"), 10)
```

```
##    seqnames     start       end  width strand     p.min  p.median
## 1     chr15  19768826  19982036 213211      * 3.073e-05 0.0013693
## 2      chr7 141419097 141441259  22163      * 9.834e-05 0.0006589
## 3     chr15  19341464  19545168 203705      * 1.407e-04 0.0012757
## 4      chr8  39356825  39497557 140733      * 2.204e-03 0.0123764
## 5      chr6  32611466  32643872  32407      * 2.220e-03 0.0150875
## 6      chr6  32059186  32065343   6158      * 2.560e-03 0.0068370
## 7     chr15  19095051  19205581 110531      * 4.307e-03 0.0150366
## 8     chr17  41785962  41914286 128325      * 6.286e-03 0.0109697
## 9      chr6  32094298  32107594  13297      * 3.235e-02 0.0462510
## 10     chr6  32066939  32093133  26195      * 3.594e-02 0.3073202
##    trans.median n.cmp
## 1        0.3304    21
## 2        0.5411     8
## 3        0.3333    30
## 4        0.5363    17
## 5        0.3614    43
## 6        0.4493     5
## 7        0.3563    13
## 8        0.4537    22
## 9        0.2576    11
## 10       0.2727    15
```

