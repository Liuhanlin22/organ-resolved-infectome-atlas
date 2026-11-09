# A Marine Fish Infectome Atlas from the South China Sea

This repository contains analysis code, figure-specific input data, and phylogenetic resources supporting the study *A Marine Fish Infectome Atlas from the South China Sea*. The analyses characterize microbial diversity, tissue-associated distributions, and host associations in wild marine fishes using metatranscriptomic data.

## Repository contents

| Directory | Contents |
| --- | --- |
| `data_script/` | R and Python scripts, input tables, and outputs organized by figure. |
| `host/` | The mitochondrial COX1 phylogenetic tree used for host identification. |
| `wildfish_phylogeny_trees/` | Microbial reference and study sequences, sequence alignments, and phylogenetic trees. |

## Requirements

Figure-generation scripts use R. Python utilities require Python 3 and the packages imported in each script. Install the R packages specified by each script before running it.

## Running the analyses

Keep the repository directory structure intact and run each plotting script from its containing directory. The corresponding input files are provided alongside the scripts or referenced through relative paths.

For example, from the repository root, generate Figure 1C with

```bash
cd data_script/figure1/figure1c
Rscript figure1C.R
```

The script saves its figure outputs in the same directory.

## Data files

Input tables are supplied in Excel, CSV, or TSV format. Sequence and alignment files are provided in FASTA format, and phylogenetic trees are provided in Newick format. Refer to the manuscript and supplementary tables for sample information and analytical methods.
