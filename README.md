# NextITS <img src='images/NextITS_logo.svg' align="right" height="70" />

![GitHub (latest release)](https://img.shields.io/github/v/release/vmikk/NextITS?label=GitHub%20release&color=23aa62)
[![Nextflow](https://img.shields.io/badge/Nextflow%20DSL2-%E2%89%A524.04-23aa62.svg)](https://www.nextflow.io/)
[![GitHub license](https://img.shields.io/github/license/vmikk/NextITS)](https://github.com/vmikk/NextITS/blob/main/LICENSE)  

[![Runs with singularity](https://img.shields.io/badge/Runs%20with-Singularity-blue?style=flat&logo=singularity)](https://sylabs.io/docs/)
[![Runs with Docker](https://img.shields.io/badge/Runs%20with-Docker-blue?style=flat&logo=docker)](https://hub.docker.com/r/vmikk/nextits/tags)  

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15074881.svg)](https://doi.org/10.5281/zenodo.15074881)

NextITS is an automated pipeline for metabarcoding fungi and other eukaryotes with full-length ITS sequenced with PacBio.
Amplicons obtained with Illumina are also supported.

## Introduction

The most widely used genetic markers for metabarcoding fungal communities are highly variable rRNA ITS1 and ITS2 sub-regions of the internal transcribed spacer. High-throughput metabarcoding has greatly improved our understanding of fungal community ecology. Here, we present NextITS, an automated pipeline for analyzing full-length ITS sequences (ITS1-5.8S-ITS2) from the Pacific Biosciences (PacBio) third-generation sequencing platform. Although the PacBio HiFi reads are highly accurate, the primary type of sequencing error is insertions or deletions in homopolymeric sites, which are also naturally common in fungal ITS. In the pipeline, we implemented correction of homopolymer errors, detection of tag-switching artefacts, and recovery of sequences false-positively annotated as chimeric. The pipeline is built using Nextflow workflow manager, with all the software dependencies packaged into Docker and Singularity containers.

## User Documentation

User documentation: https://Next-ITS.github.io/  

## Quick Start

```
nextflow run vmikk/NextITS -r main \
  -profile singularity \
  -resume \
  --input          "pacbio_ccs.fastq.gz" \
  --barcodes       "sample_barcodes.fasta" \
  --primer_forward "GTACACACCGCCCGTCG" \
  --primer_reverse "CCTSCSCTTANTDATATGC" \
  --its_region     "full" \
  --outdir         "Results"
```

## Citation

Mikryukov V., Anslan S., Tedersoo L. NextITS: a pipeline for metabarcoding fungi and other eukaryotes with full-length ITS sequenced with PacBio. [https://github.com/vmikk/NextITS](https://github.com/vmikk/NextITS)

