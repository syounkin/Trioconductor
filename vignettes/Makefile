inherited-deletions.html: inherited-deletions.asciidoc
	# echo 'library(ascii); Asciidoc("inherited-deletions.asciidoc", width=9, height=4, prefix.string="Figs/fig", driver=cacheSweaveAsciidoc(), pdf=FALSE)' | R --vanilla
	R --vanilla < asciidoc.R
	asciidoc -a data-uri -a toc -a toclevels=4 -a num inherited-deletions.asciidoc.txt
	mv inherited-deletions.asciidoc.html inherited-deletions.html

clean: 
	\rm -r *~ *.bak *.dvi *.log *.aux *.tex *.html *.asciidoc.txt
