# USUV Analysis

## Description
Data analysis of USUV sequences. Intended to generated ML phylogenetic trees for recurrent data overviews and reports. 

## Overview
This repository contains scripts and tools for analyzing Usutu virus (USUV) sequences. USUV is an arthropod-borne virus that belongs to the genus Flavivirus. This project focuses on sequence database curation, metadata curation and analysis.

## Project Structure
- Scripts for data processing and analysis
- Sequence analysis workflows

## Usage
To run the analysis scripts, ensure you have the necessary dependencies installed and execute the Shell scripts in the appropriate order.

```
cd usuv_analysis
bash main.sh --usage

```

## Getting Started
1. Clone this repository
2. Review the scripts in the project directory
3. Run the analysis scripts as needed

| Command | Description |
| :--- | :--- |
| `bash main.sh --dry-run` | Checks dependencies before starting. |
| `bash main.sh --fetch` | Downloads & curates NCBI data (incremental updates). |
| `bash main.sh --qc` | Filters out low-quality sequences (<10kb, >5% Ns). |
| `bash main.sh --merge` | Merge NCBI data with your local data. |
| `bash main.sh --alignment` | Aligns sequences using MAFFT. |
| `bash main.sh --tree` | Generates phylogeny via IQ-TREE. |
| `bash main.sh --timetree` | Run a time tree. |
| `bash main.sh --all` | Runs the entire pipeline from start to finish. |
## Contributing
Contributions are welcome! Please feel free to submit issues or pull requests.

## License
Please refer to the LICENSE file for licensing information.
