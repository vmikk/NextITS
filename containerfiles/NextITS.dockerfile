# NextITS - Dockerfile, main container
# Multi-stage build is used to compile Rust-based software
# Nextflow is included in the image

## To build the image, run:
# docker build --tag nextits --file NextITS.dockerfile .
#
## To run tests during build:
# docker build --target test --tag nextits-test --file NextITS.dockerfile .

## Build stage 1 (Rust and Cargo)
FROM rust:1.89.0-slim AS rust
RUN cargo install runiq sd

## Build stage 2 - Main
FROM rocker/r-ver:4.5.1 AS main

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV SHELL=/bin/bash
LABEL org.opencontainers.image.authors="vladimir.mikryukov@ut.ee"
LABEL org.opencontainers.image.version="1.1.0"

RUN apt-get update -qq \
  && apt-get -y --no-install-recommends install \
    tar zip unzip pigz gzip zstd xz-utils bzip2 coreutils \
    curl wget git less gawk nano rename bc \
    ca-certificates locales  \
    libtre-dev libtre5 zlib1g zlib1g-dev liblzma-dev libbz2-dev libcurl4-openssl-dev libglpk-dev libglpk40 \
    build-essential \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

## Install additional R packages
RUN install2.r --error --skipinstalled --ncpus -1 \
    remotes \
    optparse \
    R.utils \
    data.table \
    arrow \
    duckdb \
    BiocManager \
    plyr \
    dplyr \
    ggplot2 \
    doFuture \
    openxlsx \
    yaml

RUN    R -e 'BiocManager::install("Biostrings", ask = FALSE)' \
    && R -e 'BiocManager::install("ShortRead",  ask = FALSE)' \
    && R -e 'BiocManager::install("DECIPHER",   ask = FALSE)' \
    && R -e 'BiocManager::install("dada2",      ask = FALSE)' \
    && R -e 'BiocManager::install("phyloseq",   ask = FALSE)' \
    && rm -rf /tmp/downloaded_packages

RUN    install2.r --error --skipinstalled geodist phytools \
    && R -e 'ok <- tryCatch({ remotes::install_github("vmikk/metagMisc"); TRUE }, error=function(e){ message(e); FALSE }); \
          if (!ok || !requireNamespace("metagMisc", quietly=TRUE)) quit(status=1)' \
    && R -e 'ok <- tryCatch({ remotes::install_cran("qs", type = "source", configure.args = "--with-simd=AVX2"); TRUE }, error=function(e){ message(e); FALSE }); \
          if (!ok || !requireNamespace("qs", quietly=TRUE)) quit(status=1)' \
    && rm -rf /tmp/downloaded_packages

## Install conda
RUN mkdir -p /opt/software \
  && cd /opt/software \
  && curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
  && bash Miniforge3-Linux-x86_64.sh -u -b -p /opt/software/conda \
  && rm Miniforge3-Linux-x86_64.sh \
  && /opt/software/conda/bin/conda config --add channels bioconda \
  && /opt/software/conda/bin/mamba update -y --all \
  && /opt/software/conda/bin/mamba clean --all --yes

## Create conda initialization script (for Singularity compatibility)
RUN cd /opt/software \
  && { \
      echo 'eval "$(' '"/opt/software/conda/bin/conda" "shell.bash" "hook" 2> /dev/null' ')"'; \
      echo 'if [ $? -eq 0 ]; then'; \
      echo '  eval "$__conda_setup"'; \
      echo 'else'; \
      echo '  if [ -f "/opt/software/conda/etc/profile.d/conda.sh" ]; then'; \
      echo '    . "/opt/software/conda/etc/profile.d/conda.sh"'; \
      echo '  else'; \
      echo '    export PATH="/opt/software/conda/bin:$PATH"'; \
      echo '  fi'; \
      echo 'fi'; \
      echo 'unset __conda_setup'; \
    } > /opt/software/conda/init.bash

## Create conda environment and install software
RUN /opt/software/conda/bin/mamba install -y \
    "lima>=2.13.0" \
    "pbtk>=3.5.0" \
    "vsearch>=2.30.0" \
    "swarm>=3.1.5" \
    "seqkit>=2.10.1" \
    "seqfu>=1.22.3" \
    "fastp>=1.0.1" \
    "blast>=2.17.0" \
    "bioawk" \
    "miller>=6.13.0" \
    "xsv>=0.13.0" \
    "bedtools>=2.31.1" \
    "parallel>=20250622" \
    "csvtk>=0.34.0" \
    "cutadapt>=5.1" \
    "itsx>=1.1.3" \
    "bbmap>=39.33" \
    "ripgrep>=14.1.1" \
    "fd-find>=10.2.0" \
    "mmseqs2" \
  && /opt/software/conda/bin/mamba clean --all --yes


## Install cutadapt (with dependencies) from pip - it fails with conda (Python 3.13 confilict)
# RUN /opt/software/conda/bin/pip install --no-cache-dir \
#   "dnaio>=1.2.3" "xopen>=2.0.2" "cutadapt>=5.1"
  
## Add new tools (seqhasher, phredsort, ucs)
RUN cd /opt/software \
    && wget https://github.com/vmikk/seqhasher/releases/download/1.1.2/seqhasher \
    && chmod +x seqhasher \
    && mv seqhasher /opt/software/conda/bin/ \
    && wget https://github.com/vmikk/phredsort/releases/download/1.3.0/phredsort \
    && chmod +x phredsort \
    && mv phredsort /opt/software/conda/bin/ \
    && wget https://github.com/vmikk/ucs/releases/download/0.8.0/ucs \
    && chmod +x ucs \
    && mv ucs /opt/software/conda/bin/

## fqgrep
RUN git clone --depth 1 https://github.com/indraniel/fqgrep \
  && cd fqgrep \
  && make \
  && mv fqgrep /opt/software/conda/bin/ \
  && cd .. \
  && rm -r fqgrep

## rush
RUN wget https://github.com/shenwei356/rush/releases/download/v0.7.0/rush_linux_amd64.tar.gz \
  && tar -xzf rush_linux_amd64.tar.gz \
  && mv rush /opt/software/conda/bin/ \
  && rm rush_linux_amd64.tar.gz

## brename
RUN wget https://github.com/shenwei356/brename/releases/download/v2.14.0/brename_linux_amd64.tar.gz \
  && tar -xzf brename_linux_amd64.tar.gz \
  && mv brename /opt/software/conda/bin/ \
  && rm brename_linux_amd64.tar.gz

## MUMU
RUN git clone --depth 1 https://github.com/frederic-mahe/mumu.git \
    && cd ./mumu/ \
    && make && make check && make install \
  && mv mumu /opt/software/conda/bin/ \
  && cd .. \
  && rm -r mumu

## Rust tools (from the Cargo-based stage)
COPY --from=rust /usr/local/cargo/bin/runiq /opt/software/conda/bin/runiq
COPY --from=rust /usr/local/cargo/bin/sd    /opt/software/conda/bin/sd

## Update ITSx databases
RUN cd /opt/software \
    && git clone --depth 1 https://github.com/USDA-ARS-GBRU/ITS_HMMs/ \
    && find ITS_HMMs/ITSx_db/HMMs/ -name "*.hmm" | grep -v "N.hmm" \
       | /opt/software/conda/bin/parallel -j1 "/opt/software/conda/bin/hmmpress {}" \
    && rm /opt/software/conda/bin/ITSx_db/HMMs/* \
    && mv ITS_HMMs/ITSx_db/HMMs/* /opt/software/conda/bin/ITSx_db/HMMs/ \
    && rm -r ITS_HMMs \
    && sed -i '/#push(@profileSet,"Y")/s/#//' /opt/software/conda/bin/ITSx

## Install DuckDB
RUN cd /opt/software \
    && curl -L https://github.com/duckdb/duckdb/releases/download/v1.3.2/duckdb_cli-linux-amd64.zip -o duckdb_cli-linux-amd64.zip \
    && unzip duckdb_cli-linux-amd64.zip -d /opt/software/conda/bin/ \
    && rm duckdb_cli-linux-amd64.zip

## Set up environment for both Docker and Singularity compatibility
ENV PATH="/opt/software/conda/bin:${PATH}"

## Create non-privileged user  
RUN groupadd -g 1000 nextits \
    && useradd -u 1000 -g 1000 -m -s /bin/bash nextits \
    && mkdir -p /home/nextits \
    && chown -R nextits:nextits /home/nextits \
    && mkdir -p /tmp/nextits \
    && chmod 1777 /tmp/nextits

## Set software directory permissions 
## (NB! avoid recursive operations on large conda env)
RUN chmod 755 /opt/software \
    && chmod 755 /opt/software/conda \
    && chmod a+rX /opt/software/conda/bin \
    && chmod a+rX /opt/software/conda/lib

## Create entrypoint script that initializes conda properly
RUN echo '#!/bin/bash' > /opt/software/entrypoint.sh \
    && echo 'set -e' >> /opt/software/entrypoint.sh \
    && echo '# Try to source conda initialization if available' >> /opt/software/entrypoint.sh \
    && echo 'if [ -f "/opt/software/conda/init.bash" ]; then' >> /opt/software/entrypoint.sh \
    && echo '    source /opt/software/conda/init.bash' >> /opt/software/entrypoint.sh \
    && echo 'fi' >> /opt/software/entrypoint.sh \
    && echo 'exec "$@"' >> /opt/software/entrypoint.sh \
    && chmod +x /opt/software/entrypoint.sh

## Switch to non-privileged user
USER nextits
## Change working directory (for compatiblity with Singularity)
WORKDIR /tmp/nextits
ENTRYPOINT ["/opt/software/entrypoint.sh"]

## Test stage - run with: docker build --target test
FROM main AS test

# Set environment variable for R version testing
ENV R_VERSION=4.5.1

RUN echo "=== Testing R installation and packages ===" \
  && R --quiet -e "stopifnot(getRversion() == '${R_VERSION}')" \
  && echo "Testing R package installations..." \
  && printf '%s\n' \
    'required_packages <- c("optparse", "data.table", "arrow", "duckdb",' \
    '                      "plyr", "dplyr", "ggplot2", "openxlsx", "yaml",' \
    '                      "Biostrings", "DECIPHER", "dada2", "phyloseq",' \
    '                      "metagMisc", "qs")' \
    '' \
    'for(pkg in required_packages) {' \
    '  cat("Testing package:", pkg, "... ")' \
    '  tryCatch({' \
    '    suppressPackageStartupMessages(' \
    '      library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)' \
    '    )' \
    '    cat("OK\n")' \
    '  }, error = function(e) {' \
    '    cat("FAILED\n")' \
    '    stop("Package ", pkg, " failed to load: ", e$message)' \
    '  })' \
    '}' \
    'cat("All R packages loaded successfully!\n")' \
    > test_packages.R \
  && Rscript test_packages.R \
  && rm test_packages.R \
  && echo "=== Testing conda/mamba installed tools ===" \
  && tools_conda="lima bam2fastq vsearch swarm seqkit seqfu fastp blastn bioawk mlr xsv bedtools parallel csvtk ITSx cutadapt bbduk.sh rg fd mmseqs" \
  && for tool in $tools_conda; do \
       echo -n "Testing $tool... " \
       && if command -v $tool >/dev/null 2>&1; then \
            echo "OK"; \
          else \
            echo "FAILED - $tool not found in PATH" && exit 1; \
          fi; \
     done \
  && echo "=== Testing manually installed tools ===" \
  && tools_manual="seqhasher phredsort ucs fqgrep rush brename mumu duckdb runiq sd" \
  && for tool in $tools_manual; do \
       echo -n "Testing $tool... " \
       && if command -v $tool >/dev/null 2>&1; then \
            echo "OK"; \
          else \
            echo "FAILED - $tool not found in PATH" && exit 1; \
          fi; \
     done \
  && echo "=== All tests passed! Container looks ready for use ==="
