ped2fe.html: ped2fe.asciidoc
	R --vanilla < asciidoc.R
	asciidoc -a data-uri -a toc -a toclevels=4 -a num ped2fe.asciidoc.txt
	mv ped2fe.asciidoc.html ped2fe.html

clean: 
	\rm -r *~ *.bak *.dvi *.log *.aux *.tex *.html *.asciidoc.txt

