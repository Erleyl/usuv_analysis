# USUV Analysis

## Description
Data analysis of USUV sequences. Intended to generated ML phyolgenetic trees for recurrent data overviews and reports. 

## Overview
This repository contains scripts and tools for analyzing Usutu virus (USUV) sequences. USUV is an arthropod-borne virus that belongs to the genus Flavivirus. This project focuses on sequence database curation, metadata curation and analysis.

## Project Structure
- Scripts for data processing and analysis
- Sequence analysis workflows

## Usage
To run the analysis scripts, ensure you have the necessary dependencies installed and execute the Shell scripts in the appropriate order.

<cd usuv_analysis>
<bash main.sh --usage>
<bash main.sh --dry-run> to check if dependecies are properly installed
<bash main.sh --fetc> downloads all usuv sequence data available at NCBI using ncbi_datasets CLI. Next, rename sequences to a format "Accession.host.ISO_country.YYYYMMDD" and match data to a csv.file. Checks that total sequences match the metadata. The script creates a log and a "lastrun" date, so that next time is run, only data released after the prior analysis is added to the curated fasta and metadata.
<bash main.sh --qc> Runs a QC based on sequence lentgh <10kb and completeness <5% Ns. anything below will not be included
<bash main.sh --merge> Merge the NCBI data with your local data
<bash main.sh --alignment> Run MAFFT
<bash main.sh --tree> Run a phylogenetic tree using IQTREE with a Modeltest
<bash main.sh --timetree> Run a time tree
<bash main.sh --all> Runs all this steps

## Getting Started
1. Clone this repository
2. Review the scripts in the project directory
3. Run the analysis scripts as needed

## Contributing
Contributions are welcome! Please feel free to submit issues or pull requests.

## License
Please refer to the LICENSE file for licensing information.
