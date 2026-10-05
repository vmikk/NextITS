# NextITS - Dockerfile, main container
# Multi-stage build compiles Rust tools (runiq, sd) and papa2 (libpapa2.so)
# Nextflow is included in the image

## To build the image, run:
# docker build --tag nextits --file NextITS.dockerfile .
#
## To run tests during build:
# docker build --target test --tag nextits-test --file NextITS.dockerfile .

## Build stage 1.1 - Rust and Cargo
FROM rust:1.99.0-slim AS rust
RUN cargo install runiq sd

## Build stage 1.2 - papa2 (libpapa2.so)
## need to use rocker image to link against the same glibc
FROM rocker/r-ver:4.6.1 AS papa2_build
RUN apt-get update -qq \
  && apt-get -y --no-install-recommends install \
    ca-certificates git zlib1g-dev build-essential \
  && rm -rf /var/lib/apt/lists/*
WORKDIR /build
RUN git clone https://github.com/rec3141/papa2.git . \
  && git describe --tags \
  && make libpapa2.so

## Build stage 1.3 - ITSx2 source (empty by default)
## Replaced with a local directory via `--build-context itsx2_src=/path/to/itsx2`
FROM scratch AS itsx2_src

## Build stage 2 - Main
FROM rocker/r-ver:4.6.1 AS main

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV SHELL=/bin/bash
LABEL org.opencontainers.image.authors="vladimir.mikryukov@ut.ee"
LABEL org.opencontainers.image.version="1.2.4"

RUN apt-get update -qq \
  && apt-get -y --no-install-recommends install \
    tar zip unzip pigz gzip zstd xz-utils bzip2 coreutils \
    curl wget git less gawk nano rename bc \
    ca-certificates locales  \
    libtre-dev libtre5 zlib1g zlib1g-dev liblzma-dev libbz2-dev libcurl4-openssl-dev libglpk-dev libglpk40 \
    libgomp1 \
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
    yaml \
    qs2 \
    htmltools \
    jsonlite \
    vegan

RUN    R -e 'BiocManager::install("Biostrings", ask = FALSE)' \
    && R -e 'BiocManager::install("ShortRead",  ask = FALSE)' \
    && R -e 'BiocManager::install("DECIPHER",   ask = FALSE)' \
    && R -e 'BiocManager::install("dada2",      ask = FALSE)' \
    && R -e 'BiocManager::install("phyloseq",   ask = FALSE)' \
    && rm -rf /tmp/downloaded_packages

RUN install2.r --error --skipinstalled \
      geodist \
      phytools \
      ggdendro \
      gridExtra \
    && R -e 'ok <- tryCatch({ remotes::install_github("mikemc/speedyseq"); TRUE }, error=function(e){ message(e); FALSE }); \
          if (!ok || !requireNamespace("speedyseq", quietly=TRUE)) quit(status=1)' \
    && R -e 'ok <- tryCatch({ remotes::install_github("vmikk/metagMisc"); TRUE }, error=function(e){ message(e); FALSE }); \
          if (!ok || !requireNamespace("metagMisc", quietly=TRUE)) quit(status=1)' \
    && rm -rf /tmp/downloaded_packages

## Install conda
ARG MINIFORGE_VERSION=26.3.2-3
RUN mkdir -p /opt/software \
  && cd /opt/software \
  && curl -fsSL "https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-Linux-x86_64.sh" -o miniforge.sh \
  && bash miniforge.sh -b -p /opt/software/conda \
  && rm miniforge.sh \
  && /opt/software/conda/bin/conda config --system --add pinned_packages 'python=3.13' \
  && /opt/software/conda/bin/conda config --add channels bioconda \
  && /opt/software/conda/bin/mamba update -y --all \
  && /opt/software/conda/bin/mamba clean --all --force-pkgs-dirs --yes

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
    "lima>=26.2.1" \
    "pbtk>=3.5.0" \
    "vsearch>=2.32.0" \
    "swarm>=3.2.0" \
    "seqkit>=2.14.0" \
    "seqfu>=1.30.0" \
    "fastp>=1.3.7" \
    "blast>=2.17.0" \
    "bioawk" \
    "miller>=6.22.0" \
    "xsv>=0.13.0" \
    "bedtools>=2.31.1" \
    "parallel>=20260722" \
    "csvtk>=0.38.0" \
    "cutadapt>=5.2" \
    "itsx>=1.1.3" \
    "infernal>=1.1.5" \
    "bbmap>=40.02" \
    "ripgrep>=15.1.0" \
    "fd-find>=10.5.0" \
    "mmseqs2" \
  && /opt/software/conda/bin/mamba clean --all --force-pkgs-dirs --yes \
  && /opt/software/conda/bin/conda clean -afy

## papa2
RUN --mount=type=bind,from=papa2_build,source=/build,target=/tmp/papa2_src,rw \
    /opt/software/conda/bin/pip install --no-cache-dir /tmp/papa2_src \
  && SITE=$(/opt/software/conda/bin/python -c "import site; print(site.getsitepackages()[0])") \
  && ln -s papa2/libpapa2.so "$SITE/libpapa2.so" \
  && /opt/software/conda/bin/python -c "import sys, papa2, papa2.io; print('papa2 OK', papa2.__version__); \
       sys.exit(papa2.__version__.startswith('0.0.0') or not papa2.io._HAS_C_DEREP)"

  
## seqhasher, phredsort, ucs
RUN cd /opt/software \
    && wget https://github.com/vmikk/seqhasher/releases/download/1.2.0/seqhasher \
    && chmod +x seqhasher \
    && mv seqhasher /opt/software/conda/bin/ \
    && wget https://github.com/vmikk/phredsort/releases/download/1.4.0/phredsort \
    && chmod +x phredsort \
    && mv phredsort /opt/software/conda/bin/ \
    && wget https://github.com/vmikk/ucs/releases/download/0.8.0/ucs \
    && chmod +x ucs \
    && mv ucs /opt/software/conda/bin/

## USEARCH v11
RUN cd /opt/software \
    && wget https://github.com/rcedgar/usearch_old_binaries/raw/main/bin/usearch11.0.667_i86linux64 \
    && chmod +x usearch11.0.667_i86linux64 \
    && mv usearch11.0.667_i86linux64 /opt/software/conda/bin/usearch

## fqgrep
RUN git clone --depth 1 https://github.com/indraniel/fqgrep \
  && cd fqgrep \
  && make \
  && mv fqgrep /opt/software/conda/bin/ \
  && cd .. \
  && rm -r fqgrep

## rush
RUN wget https://github.com/shenwei356/rush/releases/download/v0.11.0/rush_linux_amd64.tar.gz \
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
    && git clone --depth 1 https://github.com/Mycology-Microbiology-Center/ITSx_HMMs \
    && find ITSx_HMMs/HMMs/ \( -name "*.h3f" -o -name "*.h3i" -o -name "*.h3m" -o -name "*.h3p" \) -delete \
    && find ITSx_HMMs/HMMs/ -name "*.hmm" | grep -v "N.hmm" \
       | /opt/software/conda/bin/parallel -j1 "/opt/software/conda/bin/hmmpress {}" \
    && rm /opt/software/conda/bin/ITSx_db/HMMs/* \
    && mv ITSx_HMMs/HMMs/* /opt/software/conda/bin/ITSx_db/HMMs/ \
    && rm -r ITSx_HMMs \
    && sed -i '/#push(@profileSet,"Y")/s/#//' /opt/software/conda/bin/ITSx

## ITSx2
## If a local source is supplied via `--build-context itsx2_src=...`, it is installed instead of ITSX2_SPEC
ARG ITSX2_SPEC="itsx2 @ git+https://github.com/"
RUN --mount=type=bind,from=itsx2_src,target=/tmp/itsx2_src,rw \
    if [ -f /tmp/itsx2_src/pyproject.toml ]; then \
      rm -rf /tmp/itsx2_src/build /tmp/itsx2_src/dist /tmp/itsx2_src/src/*.egg-info \
      && SPEC="/tmp/itsx2_src"; \
    else \
      SPEC="${ITSX2_SPEC}"; \
    fi \
    && /opt/software/conda/bin/pip install --no-cache-dir "${SPEC}" \
    && /opt/software/conda/bin/itsx2 itsx --version

## DuckDB
RUN cd /opt/software \
    && curl -L https://github.com/duckdb/duckdb/releases/download/v1.5.6/duckdb_cli-linux-amd64.zip -o duckdb_cli-linux-amd64.zip \
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
ENV R_VERSION=4.6.0

RUN echo "=== Testing R installation and packages ===" \
  && R --quiet -e "stopifnot(getRversion() == '${R_VERSION}')" \
  && echo "Testing R package installations..." \
  && printf '%s\n' \
    'required_packages <- c("optparse", "data.table", "arrow", "duckdb",' \
    '                      "plyr", "dplyr", "ggplot2", "openxlsx", "yaml",' \
    '                      "Biostrings", "DECIPHER", "dada2", "phyloseq",' \
    '                      "metagMisc", "qs2")' \
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
  && echo "=== Testing Python / papa2 ===" \
  && python3 -c "import papa2; print('papa2 version:', papa2.__version__)" \
  && echo "=== All tests passed! Container looks ready for use ==="
