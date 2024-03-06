# NextITS - Dockerfile, main container
# Multi-stage build is used to compile Rust-based software
# Nextflow is included in the image

## To build the image, run:
# docker build --tag nextits --file NextITS.dockerfile .

## Build stage 1 (Rust and Cargo)
FROM rust:1.76.0 AS RUST
RUN  cargo install runiq sd

## Build stage 2 - Main
FROM rocker/r-ver:4.3.2 AS MAIN

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV SHELL=/bin/bash
ENV CONDA_PREFIX="/opt/software/conda"
ENV PATH=${PATH}:"/opt/software/conda/bin/"
LABEL org.opencontainers.image.authors="vladimir.mikryukov@ut.ee"

RUN apt-get update -qq \
  && apt-get -y --no-install-recommends install \
    tar zip unzip pigz gzip xz-utils bzip2 coreutils \
    curl wget git less gawk nano rename bc \
    ca-certificates locales  \
    libtre-dev libtre5 zlib1g zlib1g-dev liblzma-dev libbz2-dev libcurl4-openssl-dev \
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
    BiocManager \
    plyr \
    ggplot2 \
    doFuture \
    openxlsx

RUN    R -e 'BiocManager::install("Biostrings", ask = FALSE)' \
    && R -e 'BiocManager::install("ShortRead", ask = FALSE)' \
    && R -e 'BiocManager::install("DECIPHER", ask = FALSE)' \
    && R -e 'BiocManager::install("dada2", ask = FALSE)' \
    && R -e 'BiocManager::install("phyloseq", ask = FALSE)' \
    && R -e 'remotes::install_github("vmikk/metagMisc")' \
    && rm -rf /tmp/downloaded_packages/

## Install conda
RUN mkdir -p /opt/software \
  && cd /opt/software \
  && curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh" \
  && bash Mambaforge-Linux-x86_64.sh -u -b -p ${CONDA_PREFIX} \
  && rm Mambaforge-Linux-x86_64.sh

## Create conda environment and install software
RUN ${CONDA_PREFIX}/bin/mamba install -y \
    -c conda-forge -c bioconda \
    "python=3.10" \
    "nextflow>=23.10.0" \
    "lima>=2.9.0" \
    "pbtk>=3.1.1" \
    "vsearch>=2.26.0" \
    "swarm>=3.1.4" \
    "seqkit>=2.6.0" \
    "seqfu>=1.20.0" \
    "fastp>=0.23.2" \
    "blast>=2.15.0" \
    "bioawk" \
    "miller>=6.10.0" \
    "bedtools>=2.31.0" \
    "parallel>=20231122" \
    "csvtk>=0.29.0" \
    "itsx>=1.1.3" \
    "cutadapt>=4.6" \
    "bbmap>=39.01" \
    "ripgrep>=14.0.0" \
    "fd-find>=8.5.3" \
    "mmseqs2" \
  && ${CONDA_PREFIX}/bin/mamba clean --all --yes

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
COPY --from=RUST /usr/local/cargo/bin/runiq /opt/software/conda/bin/runiq
COPY --from=RUST /usr/local/cargo/bin/sd    /opt/software/conda/bin/sd

WORKDIR /opt/software
