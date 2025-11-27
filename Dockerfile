FROM texlive/texlive:latest-minimal

# Add texlive bin to PATH
ENV PATH="/usr/local/texlive/2025/bin/x86_64-linux:${PATH}"

# Install fonts with Vietnamese support
RUN apt-get update && apt-get install -y --no-install-recommends \
    fonts-freefont-ttf \
    fontconfig \
    && rm -rf /var/lib/apt/lists/* \
    && fc-cache -fv

# Install xelatex and required TeX packages for Vietnamese academic documents
RUN tlmgr update --self && tlmgr install \
    xetex \
    babel-vietnamese \
    vntex \
    fontspec \
    geometry \
    setspace \
    titlesec \
    tocloft \
    fancyhdr \
    booktabs \
    multirow \
    float \
    caption \
    listings \
    xcolor \
    hyperref \
    parskip \
    enumitem \
    biblatex \
    biber \
    tools \
    latex-fonts \
    amsmath \
    amsfonts \
    ec \
    tex-gyre \
    tex-gyre-math

WORKDIR /workspace

CMD ["sh", "-c", "cd src && xelatex -interaction=nonstopmode -output-directory=../output main.tex && xelatex -interaction=nonstopmode -output-directory=../output main.tex && echo '✓ PDF generated: output/main.pdf'"]
