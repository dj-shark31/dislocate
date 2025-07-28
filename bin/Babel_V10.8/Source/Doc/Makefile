# Output directory where to write the documentation
OUTDIR=../../Doc

%.html: %.xml
	xsltproc -o $@ $<

OBJECTS=input_babel.html input_displacement.html input_patternDetect.html input_patternInit.html input_prepareDrag.html structures.html

all: 	$(OUTDIR) $(OBJECTS) INSTALL

$(OUTDIR):
	@mkdir -p $(OUTDIR)

input_babel.html:		input_xx.xsl 
input_displacement.html:	input_xx.xsl 
input_patternDetect.html:	input_xx.xsl 
input_patternInit.html:		input_xx.xsl 
input_prepareDrag.html:		input_xx.xsl 
structures.html:		structures.xsl

INSTALL:
	@rsync -au *.html $(OUTDIR)
	@rsync -au *.txt $(OUTDIR)

clean:
	@rm -f $(OBJECTS)
