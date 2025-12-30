FROM texlive/texlive:latest-small

# Install TeX packages (tex-gyre includes Times New Roman clone)
RUN tlmgr install babel-vietnamese tocloft titlesec enumitem multirow tex-gyre pgfplots biblatex biber

# Install Noto fonts for Arabic characters
RUN apt-get update && apt-get install -y fonts-noto-core && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
CMD ["sh", "-c", "cd src && xelatex -interaction=nonstopmode -halt-on-error -output-directory=../output main.tex > /dev/null && biber --output-directory=../output main && xelatex -interaction=nonstopmode -halt-on-error -output-directory=../output main.tex > /dev/null && xelatex -interaction=nonstopmode -halt-on-error -output-directory=../output main.tex | grep -E '^(!|.*:[0-9]+:|Warning|Error|Overfull|Underfull)' || true"]
