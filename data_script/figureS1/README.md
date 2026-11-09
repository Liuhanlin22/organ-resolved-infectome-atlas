# Supplementary Figure 1

| Panel | Content | Script or input |
|---|---|---|
| A | Field-identification photographs | `FigureS1A_field_identification.pdf`, three-column layout with 46 photographs |
| B | Host-genus composition within host orders | `FigureS1B_host_genus_composition.R`, `FigureS1B_host_genus_composition_input.xlsx` |
| C | Composition-partitioned log10 total-read bars above and percentage bars below | `FigureS1C_read_counts_and_percentages.R`, `FigureS1C_library_read_qc.tsv` |

## Input preparation

The source workbook is `../figure1/Table_S1_final_host_reads.xlsx`. Regenerate Figure S1B and C inputs using

```bash
python data_script/figure1/prepare_figure1_S1_inputs.py
```

Figure S1B counts each of the 134 fish individuals once. The workbook provides the proportions, individual counts, individual-level annotations and source checksum.

Figure S1C contains all 1,027 libraries. Within each host order, libraries are ordered by remaining read count from largest to smallest, with library ID breaking ties. Both plots use the same library positions. The total height of each upper bar is log10(total raw reads). Each coloured segment has height log10(total raw reads) multiplied by the original fraction of that read category. This encoding shows library size together with composition. Its definition is recorded in the plotting code. It uses the same colour definitions, category stacking order and library positions as the lower percentage bars. Each category occupies the same fraction of its library bar in both panels. The lower plot shows percentages of the total raw reads of each library. The lower x-axis labels are placed at the exact centres of their host-order blocks.

The categories are calculated as follows.

- QC-filtered reads = total raw reads minus reads retained after simple-region filtering.
- Further filtering/rRNA reads = reads retained after simple-region filtering minus actual non-rRNA input to host mapping.
- Host-mapped reads = reads in concordantly aligned pairs reported by Bowtie2.
- Remaining reads = actual non-rRNA input minus host-mapped reads.

These four categories sum to the total raw reads. The further filtering/rRNA category includes the subsequent filtering steps between simple-region filtering and the actual non-rRNA input. Its label reflects this interval. The input table also provides host-mapped and remaining percentages relative to the actual non-rRNA input.

The script provides four independent colour controls near the top: `qc_filtered_colour`, `rrna_filtered_colour`, `host_mapped_colour` and `remaining_colour`. TSV reading and writing use base R functions to support sessions with different readr versions.

## Draw the panels

```bash
Rscript data_script/figureS1/FigureS1B_host_genus_composition.R
Rscript data_script/figureS1/FigureS1C_read_counts_and_percentages.R
```

The scripts export PDF and PNG files to this directory. Figure S1C also exports `FigureS1C_read_counts_and_percentages_plot_data.tsv`, which records the exact count and percentage of every plotted category for each library. The photograph PDF is placed into the final composite during layout.

Figure S1B exports PDF with the base R PDF device and Arial font metrics so each genus name is an independent editable text object.

Figure S1C uses #AEA1E2 for further filtering/rRNA reads, #C1B77E for host-mapped reads and #CFE2BF for remaining reads. Its base-PDF export gives the legend title and each read category separate text objects.
