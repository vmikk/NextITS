# NextITS

NextITS is an automated pipeline for metabarcoding fungi and other eukaryotes with full-length ITS sequenced with PacBio.
Amplicons obtained with Illumina are also supported.

## Introduction

The most widely used genetic markers for metabarcoding fungal communities are highly variable rRNA ITS1 and ITS2 sub-regions of the internal transcribed spacer. High-throughput metabarcoding has greatly improved our understanding of fungal community ecology. Here, we present NextITS, an automated pipeline for analyzing full-length ITS sequences (ITS1-5.8S-ITS2) from the Pacific Biosciences (PacBio) third-generation sequencing platform. Although the PacBio HiFi reads are highly accurate, the primary type of sequencing error is insertions or deletions in homopolymeric sites, which are also naturally common in fungal ITS. In the pipeline, we implemented correction of homopolymer errors, detection of tag-switching artifacts, and recovery of sequences false-positively annotated as chimeric. The pipeline is built using Nextflow workflow manager, with all the software dependencies packaged into Docker and Singularity containers.
