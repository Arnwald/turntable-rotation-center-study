PYTHON ?= python3

.PHONY: all test study bem figures articles paper clean

all: study bem figures articles paper

test:
	$(PYTHON) -m pytest -q

study:
	$(PYTHON) scripts/run_study.py

bem:
	$(PYTHON) scripts/run_bem.py

figures:
	$(PYTHON) scripts/make_schematics.py
	$(PYTHON) scripts/make_figures.py
	$(PYTHON) scripts/make_bem_figures.py

articles:
	$(PYTHON) scripts/make_article_figures.py
	cd articles/audioxpress && pandoc small-turntable.md -o small-turntable.docx --resource-path=.

paper:
	cd paper && latexmk -pdf -interaction=nonstopmode main.tex

clean:
	cd paper && latexmk -C
