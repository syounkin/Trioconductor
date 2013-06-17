#!/bin/bash
name="inherited-deletions"
cd ~/trioClasses/vignettes/
pandoc -o $name.html $name.md
cp $name.html ~/sgy-website/vignettes/
rsync -a ~/trioClasses/vignettes/figures/ ~/sgy-website/vignettes/figures/
cd ~/sgy-website/
./build-website.sh
