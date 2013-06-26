inherited-deletions.html: inherited-deletions.asciidoc
	R --vanilla < asciidoc.R
	asciidoc -a data-uri -a toc -a toclevels=4 -a num inherited-deletions.asciidoc.txt
	mv inherited-deletions.asciidoc.html inherited-deletions.html

clean: 
	\rm -r *~ *.bak *.dvi *.log *.aux *.tex *.html *.asciidoc.txt
