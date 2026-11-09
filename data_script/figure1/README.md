# Figure 1

| Panel | Content | Script or input |
|---|---|---|
| A | Sampling map | Prepared during final figure assembly |
| B | COX1 host phylogeny | `figure1b/figure1B.R` and the accompanying tree, alignment and annotation files |
| C | Sequencing-library count per host species | `figure1c/figure1C.R`, `figure1c/figure1C_species_library_count.tsv` |
| D | Sampling and laboratory workflow | Prepared during final figure assembly |

`Table_S1_final_host_reads.xlsx` is the shared source for Figure 1C and Figure S1B-C. It contains 1,027 libraries, 51 host species, 134 fish individuals and 20 host orders or order-level lineages.

## Draw the panels

```bash
Rscript data_script/figure1/figure1b/figure1B.R
Rscript data_script/figure1/figure1c/figure1C.R
```

Both scripts resolve their default input and output paths relative to the script location. Figure 1C counts all libraries per host species. Its host-order blocks match the Figure S1C library ordering. The phylogenetic files were relocated from the previous Figure 1A directory to Figure 1B.
