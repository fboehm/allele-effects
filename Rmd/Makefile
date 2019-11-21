PRE1=tnseq-timbr-hotspots-all

all: $(PRE1).md $(PRE1).html 

$(PRE1).md: $(PRE1).Rmd
	R -e "rmarkdown::render('$(notdir $<)', output_format = 'github_document')"

$(PRE1).html: $(PRE1).Rmd
	R -e "rmarkdown::render('$(notdir $<)', output_format = 'html_document')"


.PHONY: all

clean:
	rm -f $(PRE1).md
	rm -f $(PRE1).html
	
