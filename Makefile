PYTHON ?= python3

.PHONY: all test study figures paper clean

all: study figures paper

test:
	$(PYTHON) -m pytest -q

study:
	$(PYTHON) scripts/run_study.py

figures:
	$(PYTHON) scripts/make_schematics.py
	$(PYTHON) scripts/make_figures.py

paper:
	cd paper && latexmk -pdf -interaction=nonstopmode main.tex

clean:
	cd paper && latexmk -C
