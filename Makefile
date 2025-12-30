IMAGE_NAME = btl-csdk-latex
OUTPUT_PDF = output/main.pdf

.PHONY: pdf build clean shell rebuild release bump-patch bump-minor bump-major

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

# Local build (if xelatex is installed) - silent mode
local:
	@cd src && xelatex -interaction=nonstopmode -halt-on-error -output-directory=../output main.tex > /dev/null
	@cd src && xelatex -interaction=nonstopmode -halt-on-error -output-directory=../output main.tex | grep -E '^(!|.*:[0-9]+:|Warning|Error|Overfull|Underfull)' || true

# Get latest tag or default to v0.0.0
LATEST_TAG := $(shell git describe --tags --abbrev=0 2>/dev/null || echo "v0.0.0")
VERSION := $(shell echo $(LATEST_TAG) | sed 's/v//')
MAJOR := $(shell echo $(VERSION) | cut -d. -f1)
MINOR := $(shell echo $(VERSION) | cut -d. -f2)
PATCH := $(shell echo $(VERSION) | cut -d. -f3)

# Version bump targets
bump-patch:
	@echo "Current version: $(LATEST_TAG)"
	@NEW_VERSION="v$(MAJOR).$(MINOR).$(shell echo $$(($(PATCH)+1)))"; \
	echo "New version: $$NEW_VERSION"; \
	git tag $$NEW_VERSION && git push origin $$NEW_VERSION

bump-minor:
	@echo "Current version: $(LATEST_TAG)"
	@NEW_VERSION="v$(MAJOR).$(shell echo $$(($(MINOR)+1))).0"; \
	echo "New version: $$NEW_VERSION"; \
	git tag $$NEW_VERSION && git push origin $$NEW_VERSION

bump-major:
	@echo "Current version: $(LATEST_TAG)"
	@NEW_VERSION="v$(shell echo $$(($(MAJOR)+1))).0.0"; \
	echo "New version: $$NEW_VERSION"; \
	git tag $$NEW_VERSION && git push origin $$NEW_VERSION

# Release (build + bump patch)
release: pdf bump-patch

# Help
help:
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  pdf        - Build PDF using Docker (default)"
	@echo "  build      - Build Docker image"
	@echo "  rebuild    - Force rebuild Docker image"
	@echo "  clean      - Remove auxiliary files"
	@echo "  distclean  - Remove all generated files"
	@echo "  shell      - Open shell in Docker container"
	@echo "  local      - Build locally (requires xelatex)"
	@echo "  bump-patch - Bump patch version (v1.0.0 -> v1.0.1) and push tag"
	@echo "  bump-minor - Bump minor version (v1.0.0 -> v1.1.0) and push tag"
	@echo "  bump-major - Bump major version (v1.0.0 -> v2.0.0) and push tag"
	@echo "  release    - Build PDF and bump patch version"
	@echo "  help       - Show this help"
