FROM texlive/texlive:latest-small

# Install TeX packages (tex-gyre includes Times New Roman clone)
RUN tlmgr install babel-vietnamese tocloft titlesec enumitem multirow tex-gyre

WORKDIR /workspace
CMD ["sh", "-c", "cd src && xelatex -interaction=nonstopmode -output-directory=../output main.tex && xelatex -interaction=nonstopmode -output-directory=../output main.tex"]
