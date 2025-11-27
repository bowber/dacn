FROM texlive/texlive:latest-small

# Install Vietnamese font and missing TeX packages
RUN apt-get update && apt-get install -y --no-install-recommends fonts-freefont-ttf \
    && rm -rf /var/lib/apt/lists/* \
    && tlmgr install babel-vietnamese tocloft titlesec enumitem multirow

WORKDIR /workspace
CMD ["sh", "-c", "cd src && xelatex -interaction=nonstopmode -output-directory=../output main.tex && xelatex -interaction=nonstopmode -output-directory=../output main.tex"]
