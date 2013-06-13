
```r
opts_chunk$set(fig.width = 7, fig.height = 5, comment = "", warning = FALSE, 
    tidy = TRUE, highlight = TRUE, fig.path = "figures/", dev = "png", background = "red")
```
















![plot of chunk hist](figures/hist.png) 





```r
c(length(DF$grange), length(reduce(DF$grange)))
```

```
[1] 632  45
```


![plot of chunk qqplot](figures/qqplot.png) 








```
   seqnames     start       end  width strand     p.min  p.median
1     chr15  19768826  19982036 213211      * 3.073e-05 0.0013693
2      chr7 141419097 141441259  22163      * 9.834e-05 0.0006589
3     chr15  19341464  19545168 203705      * 1.407e-04 0.0012757
4      chr8  39356825  39497557 140733      * 2.204e-03 0.0123764
5      chr6  32611466  32643872  32407      * 2.220e-03 0.0150875
6      chr6  32059186  32065343   6158      * 2.560e-03 0.0068370
7     chr15  19095051  19205581 110531      * 4.307e-03 0.0150366
8     chr17  41785962  41914286 128325      * 6.286e-03 0.0109697
9      chr6  32094298  32107594  13297      * 3.235e-02 0.0462510
10     chr6  32066939  32093133  26195      * 3.594e-02 0.3073202
11     chr6  32650822  32664356  13535      * 3.960e-02 0.0660980
12    chr16  33778130  33820307  42178      * 4.223e-02 0.3164739
13    chr11  55124465  55209499  85035      * 5.743e-02 0.5299039
14     chr3  75502426  75719139 216714      * 5.824e-02 0.3020938
15    chr14  18347035  18372086  25052      * 5.877e-02 0.4131213
16     chr5  97074222  97125076  50855      * 7.163e-02 0.0854169
17     chr1 103941535 104099390 157856      * 1.051e-01 0.2673147
18     chr5  69359352  69433008  73657      * 1.058e-01 0.1806520
19     chr6  31388080  31397263   9184      * 1.339e-01 0.3760968
20    chr19  20404485  20507068 102584      * 1.457e-01 0.3730939
21    chr11  48890168  48918267  28100      * 1.614e-01 0.4851230
22    chr16  32404517  32530051 125535      * 1.907e-01 0.6826270
23    chr12  36404411  36532019 127609      * 2.296e-01 0.3155882
24     chr1 195087039 195087039      1      * 2.584e-01 0.2583597
25     chr9  43594114  43674189  80076      * 2.810e-01 0.4630094
   trans.median n.cmp
1        0.3304    21
2        0.5411     8
3        0.3333    30
4        0.5363    17
5        0.3614    43
6        0.4493     5
7        0.3563    13
8        0.4537    22
9        0.2576    11
10       0.2727    15
11       0.2239     9
12       0.5162    12
13       0.5346    26
14       0.4319    15
15       0.4548    12
16       0.4783     5
17       0.4800    25
18       0.2623     9
19       0.4776    14
20       0.5013     6
21       0.4653    20
22       0.2769   123
23       0.1846     7
24       0.4035     1
25       0.3715    14
```


![plot of chunk transvp](figures/transvp.png) 



```
GRanges with 1 range and 4 metadata columns:
      seqnames               ranges strand |             p.min
         <Rle>            <IRanges>  <Rle> |         <numeric>
  [1]    chr16 [32404517, 32530051]      * | 0.190686032146876
               p.median      trans.median     n.cmp
              <numeric>         <numeric> <integer>
  [1] 0.682627027049546 0.276872964169381       123
  ---
  seqlengths:
            chr1   chr1_random          chr2 ...          chrY          chrM
       247249719       1663265     242951149 ...      57772954         16571
```

The outlier is on chromsome 16.  It is a region with \Sexpr{values(bad.region.gr)$n.cmp} components, and has width \Sexpr{width(bad.region.gr)/1e3} kB. chr16:\Sexpr{start(bad.region.gr)}-\Sexpr{end(bad.region.gr)}.  If we remove the outlying region on chromsome 16 we see the following.


![plot of chunk transvp2](figures/transvp2.png) 



![plot of chunk cumsum](figures/cumsum.png) 


![plot of chunk transmedianhist](figures/transmedianhist.png) 


```r
thresh <- with(as(values(regions.gr), "data.frame"), median(trans.median))
regions.gr.clean <- regions.gr[which(values(regions.gr)$trans.median >= thresh)]
DF.clean <- DF[queryHits(findOverlaps(DF$grange, regions.gr.clean)), ]
```


![plot of chunk qqplot-clean](figures/qqplot-clean.png) 



![plot of chunk cumsum2](figures/cumsum2.png) 



![plot of chunk phist](figures/phist.png) 



![plot of chunk transvp3](figures/transvp3.png) 
