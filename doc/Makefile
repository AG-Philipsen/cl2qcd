## LaTeX Makefile

MAINDOC=main

TEX=latex
PDFTEX=pdflatex
BIB=bibtex

DVIPS=dvips

TEXSRC = $(wildcard *.tex)

FIGSRC = $(wildcard figures/*.eps) \
         $(wildcard figures/*.jpg) \

all: ps

ps: $(MAINDOC).ps
pdf: $(MAINDOC).pdf

gv: ps
	gv $(MAINDOC).ps
kpdf: pdf
	kpdf $(MAINDOC).pdf
evince: pdf
	evince $(MAINDOC).pdf

clean:
	rm -f *.aux *.bbl *.blg *.log *.toc

cleancomplete: clean
	rm -f $(MAINDOC).{dvi,ps,pdf}

$(MAINDOC).tex: $(TEXSRC) $(FIGSRC)
	$(TEX) $(MAINDOC)
$(MAINDOC).aux: $(TEXSRC) $(FIGSRC)
	$(TEX) $(MAINDOC)
$(MAINDOC).bbl: $(TEXSRC) $(FIGSRC) $(BIBFILE)
	$(TEX) $(MAINDOC)
	$(BIB) $(MAINDOC)
	$(TEX) $(MAINDOC)
$(MAINDOC).dvi: $(TEXSRC) $(FIGSRC) $(MAINDOC).bbl $(MAINDOC).aux
	$(TEX) $(MAINDOC)
$(MAINDOC).ps: $(MAINDOC).dvi
	$(DVIPS) $(MAINDOC)
$(MAINDOC).pdf: $(TEXSRC) $(FIGSRC) $(MAINDOC).bbl $(MAINDOC).aux
	$(PDFTEX) $(MAINDOC)

