#!/bin/bash
# name="inherited-deletions"

name="rawData"
cd ~/trioClasses/vignettes/
pandoc -o $name.html $name.md
cp $name.html ~/sgy-website/vignettes/
rsync -a ~/trioClasses/vignettes/figures/ ~/sgy-website/vignettes/figures/

name="rawData-chr7"
cd ~/trioClasses/vignettes/
pandoc -o $name.html $name.md
cp $name.html ~/sgy-website/vignettes/
rsync -a ~/trioClasses/vignettes/figures/ ~/sgy-website/vignettes/figures/

name="rawData-chr15"
cd ~/trioClasses/vignettes/
pandoc -o $name.html $name.md
cp $name.html ~/sgy-website/vignettes/
rsync -a ~/trioClasses/vignettes/figures/ ~/sgy-website/vignettes/figures/

name="rawData-chr6"
cd ~/trioClasses/vignettes/
pandoc -o $name.html $name.md
cp $name.html ~/sgy-website/vignettes/
rsync -a ~/trioClasses/vignettes/figures/ ~/sgy-website/vignettes/figures/

cd ~/sgy-website/
./build-website.sh

