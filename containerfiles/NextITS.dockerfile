# NextITS - Dockerfile, main container
# Multi-stage build is used to compile Rust-based software
# Nextflow is included in the image

## To build the image, run:
# docker build --tag nextits --file NextITS.dockerfile .

## Build stage 1 (Rust and Cargo)
FROM rust:1.84.1 AS rust
RUN cargo install runiq sd

## Build stage 2 - Main
FROM rocker/r-ver:4.4.2 AS main

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV SHELL=/bin/bash
ENV CONDA_PREFIX="/opt/software/conda"
ENV PATH=${PATH}:"/opt/software/conda/bin/"
LABEL org.opencontainers.image.authors="vladimir.mikryukov@ut.ee"
LABEL org.opencontainers.image.version="0.8.2"

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
    openxlsx

RUN    R -e 'BiocManager::install("Biostrings", ask = FALSE)' \
    && R -e 'BiocManager::install("ShortRead", ask = FALSE)' \
    && R -e 'BiocManager::install("DECIPHER", ask = FALSE)' \
    && R -e 'BiocManager::install("dada2", ask = FALSE)' \
    && R -e 'BiocManager::install("phyloseq", ask = FALSE)' \
    && R -e 'remotes::install_github("vmikk/metagMisc")' \
    && R -e 'remotes::install_cran("qs", type = "source", configure.args = "--with-simd=AVX2")' \
    && rm -rf /tmp/downloaded_packages/

## Install conda
RUN mkdir -p /opt/software \
  && cd /opt/software \
  && curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  && bash Miniconda3-latest-Linux-x86_64.sh -u -b -p ${CONDA_PREFIX} \
  && rm Miniconda3-latest-Linux-x86_64.sh \
  && ${CONDA_PREFIX}/bin/conda config --add channels defaults \
  && ${CONDA_PREFIX}/bin/conda config --add channels conda-forge \
  && ${CONDA_PREFIX}/bin/conda config --add channels bioconda \
  && ${CONDA_PREFIX}/bin/conda install --quiet --yes mamba "python=3.12"

## Create conda environment and install software
RUN ${CONDA_PREFIX}/bin/mamba install -y \
    -c conda-forge -c bioconda \
    "lima>=2.12.0" \
    "pbtk>=3.4.0" \
    "vsearch>=2.29.4" \
    "swarm>=3.1.5" \
    "seqkit>=2.9.0" \
    "seqfu>=1.22.3" \
    "fastp>=0.24.0" \
    "blast>=2.16.0" \
    "bioawk" \
    "miller>=6.13.0" \
    "xsv>=0.13.0" \
    "bedtools>=2.31.1" \
    "parallel>=20241122" \
    "csvtk>=0.31.0" \
    "itsx>=1.1.3" \
    "cutadapt>=5.0" \
    "bbmap>=39.15" \
    "ripgrep>=14.1.1" \
    "fd-find>=10.2.0" \
    "mmseqs2" \
  && ${CONDA_PREFIX}/bin/mamba clean --all --yes

## Add new tools (seqhasher, phredsort, ucs)
RUN cd /opt/software \
    && wget https://github.com/vmikk/seqhasher/releases/download/1.1.1/seqhasher \
    && chmod +x seqhasher \
    && mv seqhasher ${CONDA_PREFIX}/bin/ \
    && wget https://github.com/vmikk/phredsort/releases/download/1.3.0/phredsort \
    && chmod +x phredsort \
    && mv phredsort ${CONDA_PREFIX}/bin/ \
    && wget https://github.com/vmikk/ucs/releases/download/0.8.0/ucs \
    && chmod +x ucs \
    && mv ucs ${CONDA_PREFIX}/bin/

## fqgrep
RUN git clone --depth 1 https://github.com/indraniel/fqgrep \
  && cd fqgrep \
  && make \
  && mv fqgrep ${CONDA_PREFIX}/bin/ \
  && cd .. \
  && rm -r fqgrep

## rush
RUN wget https://github.com/shenwei356/rush/releases/download/v0.5.4/rush_linux_amd64.tar.gz \
  && tar -xzf rush_linux_amd64.tar.gz \
  && mv rush ${CONDA_PREFIX}/bin/ \
  && rm rush_linux_amd64.tar.gz

## brename
RUN wget https://github.com/shenwei356/brename/releases/download/v2.14.0/brename_linux_amd64.tar.gz \
  && tar -xzf brename_linux_amd64.tar.gz \
  && mv brename ${conda_prefix}/bin/ \
  && rm brename_linux_amd64.tar.gz

## MUMU
RUN git clone --depth 1 https://github.com/frederic-mahe/mumu.git \
    && cd ./mumu/ \
    && make && make check && make install \
  && mv mumu ${CONDA_PREFIX}/bin/ \
  && cd .. \
  && rm -r mumu

## Rust tools (from the Cargo-based stage)
COPY --from=rust /usr/local/cargo/bin/runiq /opt/software/conda/bin/runiq
COPY --from=rust /usr/local/cargo/bin/sd    /opt/software/conda/bin/sd

## Update ITSx databases
RUN cd /opt/software \
    && git clone --depth 1 https://github.com/USDA-ARS-GBRU/ITS_HMMs/ \
    && find ITS_HMMs/ITSx_db/HMMs/ -name "*.hmm" | grep -v "N.hmm" \
       | ${CONDA_PREFIX}/bin/parallel -j1 "${CONDA_PREFIX}/bin/hmmpress {}" \
    && rm ${CONDA_PREFIX}/bin/ITSx_db/HMMs/* \
    && mv ITS_HMMs/ITSx_db/HMMs/* ${CONDA_PREFIX}/bin/ITSx_db/HMMs/ \
    && rm -r ITS_HMMs \
    && sed -i '/#push(@profileSet,"Y")/s/#//' ${CONDA_PREFIX}/bin/ITSx

## Install DuckDB
RUN cd /opt/software \
    && curl -L https://github.com/duckdb/duckdb/releases/download/v1.2.0/duckdb_cli-linux-amd64.zip -o duckdb_cli-linux-amd64.zip \
    && unzip duckdb_cli-linux-amd64.zip -d ${CONDA_PREFIX}/bin/ \
    && rm duckdb_cli-linux-amd64.zip

WORKDIR /opt/software
