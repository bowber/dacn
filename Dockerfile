FROM texlive/texlive:latest-small

# Install TeX packages (tex-gyre includes Times New Roman clone)
RUN tlmgr install babel-vietnamese tocloft titlesec enumitem multirow tex-gyre

# Install Noto fonts for Arabic characters
RUN apt-get update && apt-get install -y fonts-noto-core && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
CMD ["sh", "-c", "cd src && xelatex -interaction=nonstopmode -output-directory=../output main.tex && xelatex -interaction=nonstopmode -output-directory=../output main.tex"]
