# AGENTS.md - AI Coding Agent Instructions

## Build Commands
- `make pdf` - Build PDF using Docker (default, recommended)
- `make local` - Build locally (requires xelatex)
- `make clean` - Remove auxiliary files (.aux, .log, .toc, etc.)
- `make distclean` - Remove all generated files including PDFs
- `make shell` - Open interactive shell in Docker container for debugging

## Versioning
- Versioned via git tags (e.g., `v1.0.11`)
- Use semantic versioning: `vMAJOR.MINOR.PATCH`
- Only create new version tags when explicitly requested by the user
- Commit changes after completing each task/feature

## Project Structure
LaTeX academic lab report using XeLaTeX with Vietnamese language support.
- `src/` - LaTeX source files (main.tex is entry point, includes other .tex files)
- `assets/` - Images and logos
- `data/` - Experimental data (.csv, .TRD files)
- `tmp/` - Temporary files (git-ignored, use for extracted images, intermediate files, etc.)

## Code Style Guidelines
- **Section headers**: Use `% ============================================================` comment blocks
- **Subsection headers**: Use `% ------------------------------------------------------------` comment blocks
- **Labels**: Use snake_case with descriptive prefixes (e.g., `fig:feedback_control`)
- **Language**: Vietnamese text throughout; use `babel-vietnamese` package features
- **TikZ diagrams**: Use predefined styles in preamble.tex (`block`, `sum`, `arrow`, `line`)
- **Modular structure**: Keep each lab in separate file (lab2.tex, lab3.tex, etc.)
- **Adding/removing labs**: Update `\input{}` in main.tex; ToC auto-generates from sections
- **Table captions**: Place `\caption{}` AFTER `\end{tabular}` (below the table)
