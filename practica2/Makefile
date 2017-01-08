PDF=$(patsubst %.md, %.pdf, $(wildcard *.md))

default: $(PDF)

%.pdf: %.md
	pandoc $< -o $@
