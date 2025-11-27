FROM texlive/texlive:latest-small

# Install Times New Roman (Microsoft fonts) and missing TeX packages
RUN echo "deb http://deb.debian.org/debian testing contrib" >> /etc/apt/sources.list \
    && echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections \
    && apt-get update && apt-get install -y --no-install-recommends ttf-mscorefonts-installer \
    && rm -rf /var/lib/apt/lists/* \
    && fc-cache -f \
    && tlmgr install babel-vietnamese tocloft titlesec enumitem multirow

WORKDIR /workspace
CMD ["sh", "-c", "cd src && xelatex -interaction=nonstopmode -output-directory=../output main.tex && xelatex -interaction=nonstopmode -output-directory=../output main.tex"]
