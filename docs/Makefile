# make | grep -A 3 -iE '!|Warning|Error|failed|Success|dvips:'
OUT=transforms
SRC=$(OUT).tex

LATEX=latex
LATEXFLAGS=-interaction=nonstopmode
DVIPS=dvips
DVIPSFLAGS=
PDFLATEX=pdflatex
PDFLATEXFLAGS=-synctex=1 -interaction=nonstopmode
BIBTEX=bibtex
PS2PDF=ps2pdf
SCCSS=Success

.PHONE: all latex pdflatex bibtex ps2pdf clean

all: latex

latex:
	@$(LATEX) $(LATEXFLAGS) $(SRC) && $(DVIPS) -q -o $(OUT).ps $(OUT).dvi && echo Success

pdflatex:
	$(PDFLATEX) $(PDFLATEXFLAGS) $(SRC)

bibtex:
	$(BIBTEX) $(OUT).aux

ps2pdf:
	$(PS2PDF) $(OUT).ps

clean:
	rm -f $(OUT).log $(OUT).pdf $(OUT).dvi $(OUT).ps $(OUT).aux $(OUT).out $(OUT).synctex.gz
