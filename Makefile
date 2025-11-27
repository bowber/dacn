IMAGE_NAME = btl-csdk-latex
OUTPUT_PDF = output/main.pdf

.PHONY: pdf build clean shell rebuild

# Build PDF using Docker
pdf: build
	docker run --rm -v $(PWD):/workspace $(IMAGE_NAME)

# Build Docker image
build:
	@docker build -t $(IMAGE_NAME) . 2>/dev/null || docker build -t $(IMAGE_NAME) .

# Force rebuild Docker image
rebuild:
	docker build --no-cache -t $(IMAGE_NAME) .

# Clean auxiliary files
clean:
	rm -f output/*.aux output/*.log output/*.toc output/*.out output/*.fls output/*.fdb_latexmk

# Deep clean (including PDF)
distclean: clean
	rm -f output/*.pdf

# Interactive shell in container
shell:
	docker run --rm -it -v $(PWD):/workspace $(IMAGE_NAME) /bin/bash

# Local build (if xelatex is installed)
local:
	cd src && xelatex -interaction=nonstopmode -output-directory=../output main.tex
	cd src && xelatex -interaction=nonstopmode -output-directory=../output main.tex

# Help
help:
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  pdf      - Build PDF using Docker (default)"
	@echo "  build    - Build Docker image"
	@echo "  rebuild  - Force rebuild Docker image"
	@echo "  clean    - Remove auxiliary files"
	@echo "  distclean- Remove all generated files"
	@echo "  shell    - Open shell in Docker container"
	@echo "  local    - Build locally (requires xelatex)"
	@echo "  help     - Show this help"
