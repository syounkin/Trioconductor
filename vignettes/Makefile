Trioconductor.html: Trioconductor.asciidoc
	R --vanilla < asciidoc.R
	asciidoc -a data-uri -a toc -a toclevels=4 -a num Trioconductor.asciidoc.txt
	mv Trioconductor.asciidoc.html Trioconductor.html
	rm Trioconductor.asciidoc.txt

clean: 
	\rm -r *~ *.bak *.dvi *.log *.aux *.tex *.html *.asciidoc.txt

