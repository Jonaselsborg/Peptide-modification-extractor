# Quality Filtering of ADP-ribosylation Site Data from MaxQuant

This repository provides an R script for processing **MaxQuant evidence.txt** files containing ADP-ribosylation (ADPr) site information.  
The script extracts localization scores, applies filtering to retain only high-confidence sites, annotates each site within the corresponding UniProt ID, and transfers match-between-runs (MBR) intensities for modified peptides if they were detected by MS/MS.

## Requirements

- **File naming convention**: Replicates must have the suffix `_01` (underscore followed by two digits), while the rest of the filename should remain identical. This convention is used to group replicates into experimental sets.  
- **R packages**:
  - [`here`](https://CRAN.R-project.org/package=here) – for defining the working directory.  
  - [`seqinr`](https://CRAN.R-project.org/package=seqinr) – for reading FASTA files.  
  - [`tidyverse`](https://www.tidyverse.org/) – for data wrangling and general utilities.  

Ensure these packages are installed prior to running the script.

## Input Files

Two files are required in the root of the repository:  

1. **evidence.txt** – The MaxQuant evidence file containing ADPr site identifications.  
2. **FASTA file** – The reference FASTA file used during the MaxQuant database search (extension: `.fasta`).  

> ⚠️ Note: Ensure that the FASTA file corresponds exactly to the one used in your MaxQuant search.

## Script Configuration

The script currently expects two specific columns from the evidence file:  

- `ADP-ribosylation (CDEHKRSTY)`  
- `ADP-ribosylation (CDEHKRSTY) Probabilities`  

Depending on how your MaxQuant search is configured, these column names may differ. If so, please rename them to match the expected format before running the script.

## Running the Script

You can run the script in two ways:  

- **From RStudio** – Open the script file and execute it interactively.  
- **From the terminal** – Run the script with:  
  ```bash
  Rscript ADPr-evidence_extractor.R
