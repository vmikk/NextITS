# Changelog

All notable changes to this project will be documented in this file.  

This project tries to adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).  
For version numbering, we use the following convention: `MAJOR.MINOR.PATCH`.  
Each element increases numerically (e.g., `1.9.0` -> `1.10.0` -> `1.11.0`).  


## [1.0.0] - 2025-03-24

- Added support of asymmetric barcoding scheme for demultiplexing of PacBio data  
- Added support of BAM files (CCS) as input  
- Added support for SWARM *d*=1 pre-clustering  
- Changed the selection of representative sequences (sequence with the highest quality score is taken as the representative; using [`phredsort`](https://github.com/vmikk/phredsort))  
- Refactored sequence quality estimation  
- Improved processing speed (using DuckDB and Parquet format)  
- Improved tag valiadion for demultiplexing  
- Improved compression speed for output files (runs in parallel using `pigz`)  
- Update of the database for reference-based chimera detection (using the [EUKARYOME database](https://eukaryome.org/))  
- New parameters added:  
    - `step` (specifies which pipeline step to run - "Step1" or "Step2")  
    - `storagemode` (Adjusts how files are directed to the results folder)  
    - `gzip_compression` (Controls GZIP compression level in output files)  
    - categorical `lima_barcodetype` replaces boolean `lima_dualbarcode`  
    - `lima_minendscore` (For asymmetric and dual barcoding scheme)  
    - `lima_minrefspan` (Controls barcode coverage)  
    - `lima_minscoringregions` (Controls the number of reqired barcodes for demultiplexing using dual barcodes)  
- Added auxilarry output files:  
    - [Step-1] All rRNA parts extracted by ITSx (pooled within sequencing run - useful for extracting these regions for representative sequences)  
    - [Step-1] File with quality scores for full-length sequences (after QC and trimming)  
    - [Step-2] File with joined sequence memebership (dereplication, pre-clustering, and clustering)  
- Primer trimming prior ITSx is now default (sequence quality is also estimated on trimmed sequence)  
- Fixed VSEARCH clustering on denoised reads  
- Resolved an issue where no *de novo* chimeras were detected  
- Reconfigured parameter specification  
- Introduced a parameter schema and enhanced parameter validation  
- Container updates to included the latest versions of dependencies  
- New dependencies - specialized tools written in Go to speed up the processing:  
    - [`phredsort` (https://github.com/vmikk/phredsort)](https://github.com/vmikk/phredsort) (Sorts sequences by quality score)  
    - [`seqhasher` (https://github.com/vmikk/seqhasher)](https://github.com/vmikk/seqhasher) (Hashes sequences)  
    - [`ucs` (https://github.com/vmikk/ucs)](https://github.com/vmikk/ucs) (Parses UC files and converts them to parquet format)  


## [0.5.0] - 2023-08-08

- New `seqstats` sub-workflow (only dereplication, primer validation, and basic run stats)  
- Add SWARM clustering ([Mahé et al., 2022 DOI:10.1093/bioinformatics/btab493](https://academic.oup.com/bioinformatics/article/38/1/267/6318385))  
- Add post-clustering curation with LULU ([Frøslev et al., 2017 DOI:doi.org/10.1038](https://www.nature.com/articles/s41467-017-01312-x))  
- Add barcode validation step  
- Add SSU and LSU region-based output sequences  
- Add support for UNOISE-only output (without clustering)  
- Add `merge_replicates` parameter (Step-2) for merging or keeping separate sample replicates  
- Update Step-1 run summary (add homopolymer stats)  
- Deprecate taxonomy annotation workflow at Step-1  
- Fixed different extensions in demultiplexed input  
- Experimental: UNITE-style dereplication (allows query sequences to vary in length at 100% similarity)  
- Experimental: support of alternative alignment penalty scores (ITS-specific feature)


## [0.4.0] - 2023-05-08

- Add Step-2 workflow for pooling, dereplicating, and clustering sequences from Step-1  
    - Read clustering with VSEARCH ([Rognes et al., 2016 DOI:10.7717/peerj.2584](https://peerj.com/articles/2584/))  
    - Error-correction with UNOISE2 ([Edgar, 2016 DOI:10.1101/081257](https://www.biorxiv.org/content/10.1101/081257v1))  
- Add run summary for Step-1 (read counts at different pipeline stages)  
- Separate config for HPC clusters  
- Add Docker container  


## [0.3.0] - 2023-03-02

- Add support for pre-demultiplexed data as input  
- Add option for semi-full-length ITS (especially useful when forward primer is located at the very end of SSU and the HMM site can not be recognized by ITSx)  
- Add removal of long homopolymer artefacts at QC stage  
- Correct handling of a case with no valid sequences at primer checking step (thank to Taavi Riit for reporting the bug)  
- Bug fixed in `assemble_its` (thanks to Kadri Põldmaa for discovering the error)  
- Addition of ITSx detailed results (with information on the HMM profile used for ITS extraction)  
- Fixed sample names for the rescued chimeric sequences  
- Minor fixes related with the Singularity container, output directory, help message, and single-end QC  

## [0.2.0] - 2022-09-30

- Add Ilumina-based workflow (see `--seqplatform` flag)  
- Publish Singularity image to Singularity library  
- Minor bugfixes in `primer_check` (multiprimer artefacts), `pool_seqs` (sequence headers), and `prep_asvtab` (aggregation of non-unique joined Illumina sequences) processes  
- New logo design (thanks to Olesya Dulya)  

## [0.0.1] - 2022-07-07

- Initial release  
