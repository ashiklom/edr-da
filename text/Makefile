.PHONY: all clean

all:
	latexmk -pdf main.tex
	latexpand main.tex > main-expanded.tex

diff:
	latexdiff --config="PICTUREENV=(?:picture|DIFnomarkup|align)[\w\d*@]*" \
		--math-markup=0 \
		--disable-citation-markup \
		--flatten \
		../_preprint-version/text/main.tex main.tex > diff.tex
	latexmk -pdf diff.tex

clean:
	latexmk -c
