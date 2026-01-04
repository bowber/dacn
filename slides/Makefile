# Makefile for Neubrutalism Beamer Slides
# Uses Docker for consistent builds

.PHONY: all clean

# Default target
all: main.pdf

# Build PDF using Docker
main.pdf: main.tex preamble.tex
	cd .. && docker run --rm -v $(PWD)/..:/workspace -w /workspace/slides btl-csdk-latex xelatex -output-directory=. main.tex
	cd .. && docker run --rm -v $(PWD)/..:/workspace -w /workspace/slides btl-csdk-latex xelatex -output-directory=. main.tex

# Clean auxiliary files
clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb

# Clean all including PDF
distclean: clean
	rm -f main.pdf
