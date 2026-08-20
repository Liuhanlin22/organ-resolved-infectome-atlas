# Wild Fish Pathogen Analysis, related to An Organ-Resolved Infectome Atlas of Wild Fishes

This repository contains the analysis scripts, input tables, phylogenetic files and figure outputs for the wild aquatic animal pathogen study. The scripts are intended for reproducible figure generation and use paths relative to their own directories.

## Repository layout

```text
github_upload/
├── data_script/
│   ├── figure1/          Host phylogeny and library quality control
│   ├── figure2/          Pathogen overview, heatmaps and abundance summaries
│   ├── figure3/          Organ-specific pathogen abundance plots
│   ├── figure4/          Organ distribution, prevalence, richness and PCoA
│   ├── figure5/          Host-level summaries, ordination and co-infection analyses
│   ├── figure6/          Host-pathogen network and updated heatmap analyses
│   ├── figureS1/         Viral host-genus composition
│   ├── figureS4/         Organ preference and complete-eight-organ analyses
│   ├── figureS5/         Pathogen-by-individual UpSet analysis
│   └── figureS6/         Updated host-species heatmap and cross-range summaries
├── host/                 Host phylogenetic input files
└── wildfish_phylogeny_trees/
    ├── Broad-scale trees, related to Figure S2/
    └── Fine-scale trees, related to Figure S3/
```

Each figure directory contains the relevant R script, input files and, where available, the generated PDF, SVG, PNG or TIFF figures and machine-readable result tables.

## Analysis conventions

- The primary abundance matrix contains 75 pathogen taxa and 1,027 libraries.
- Organ-level paired analyses use 112 fish with all eight organs successfully sampled and sequenced, corresponding to 896 libraries.
- Complete-organ analyses treat fish identity as the independent biological unit and organs as repeated measurements.

## Software requirements

The analyses were prepared for R 4.4.2. Install the packages required by the figure being reproduced. The complete package set used across the scripts is:

```r
install.packages(c(
  "openxlsx", "readxl", "writexl", "ggplot2", "dplyr", "tidyr",
  "readr", "tibble", "patchwork", "scales", "gridExtra", "pheatmap",
  "vegan", "permute", "Rtsne", "VennDiagram", "igraph", "ggraph",
  "ggrepel", "cowplot", "RColorBrewer", "networkD3", "htmlwidgets",
  "webshot", "svglite", "ragg", "jsonlite", "ape"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ggtree")
```

The Sankey workflow may additionally require a local PhantomJS installation for `webshot`. Some scripts use optional packages only for specific output formats or font/rendering devices.

The helper script `data_script/figure6/figure6A/prepare_Figure6A_input.py` requires Python 3 with `pandas`, `numpy` and `openpyxl`.

## Running a figure script

Run a script from the repository root or from its own figure directory. For example:

```bash
Rscript data_script/figure4/Figure4E_richness_complete8.R
Rscript data_script/figureS4/figureS4B_C/FigureS4B_prepare_and_plot.R
Rscript data_script/figureS4/figureS4B_C/figureS4C/FigureS4C_pathogen_abundance_by_tissue_paired.R
Rscript data_script/figureS6/figureS6A/FigureS6A_updated.R
Rscript data_script/figureS6/FigureS6B_C_updated.R
```

Scripts read the input files stored with the corresponding figure directory and write formal result tables and figure files to that directory. Existing figure files may be overwritten when a script is rerun.

## Input and output formats

- `.xlsx`: abundance matrices, sample metadata, pathogen annotations and statistical result tables.
- `.csv`, `.tsv`: plotting data, annotations, summaries, network nodes and edges.
- `.R`: figure-generation and statistical-analysis scripts.
- `.py`: optional input-preparation helper for Figure 6A.
- `.pdf`, `.svg`: publication-quality vector figures where provided.
- `.png`, `.tiff`: raster previews or high-resolution figure exports.
- `.nwk`, `.fas`: phylogenetic trees, alignments and sequence data used for Figures 1, S2 and S3.

## Reproducibility notes

Run scripts with the input files supplied in this repository. Do not replace the complete-eight-organ matrix with the 1,027-library matrix for analyses explicitly labeled as complete-organ or paired analyses. Conversely, Figure S6 analyses use their supplied updated matrix and annotation workbook. Species names should be kept consistent with the current Table S1 and Figure 6A annotation.

The repository contains generated figure files for direct inspection. The R scripts and machine-readable tables are the authoritative sources for reproducing the analyses.

## License and data use

The analysis code in this repository is released under the MIT License. 
