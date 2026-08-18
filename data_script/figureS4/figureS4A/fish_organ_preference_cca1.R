#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pheatmap)
  library(vegan)
  library(RColorBrewer)
})

working_candidate <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_name <- "fish_individual_organ_rpm.csv"
script_dir <- if (file.exists(file.path(working_candidate, input_name))) {
  working_candidate
} else {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) != 1) stop("Run this script from its own directory.")
  dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/"))
}

rpm_file <- file.path(script_dir, input_name)
metadata_file <- file.path(script_dir, "fish_lib_info_cca.csv")
taxonomy_file <- file.path(script_dir, "fish_taxon.csv")
output_pdf <- file.path(script_dir, "FigureS4A_CCA_heatmap.pdf")
output_png <- file.path(script_dir, "FigureS4A_CCA_heatmap_600dpi.png")
setwd(script_dir)
output_pdf_device <- basename(output_pdf)
output_png_device <- basename(output_png)

rpm <- read.csv(rpm_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
taxonomy <- read.csv(taxonomy_file, row.names = 1, check.names = FALSE)

if (!identical(rownames(rpm), rownames(metadata))) {
  stop("Sample IDs differ between RPM and metadata files.")
}
if (!identical(colnames(rpm), rownames(taxonomy))) {
  stop("Pathogen names/order differ between RPM and taxonomy files.")
}
if (!all(vapply(rpm, is.numeric, logical(1)))) stop("All RPM columns must be numeric.")
if (any(!is.finite(as.matrix(rpm))) || any(as.matrix(rpm) < 0)) {
  stop("RPM values must be finite and non-negative.")
}
if (!identical(names(metadata), "organ")) stop("Metadata must contain one column named organ.")

names(metadata) <- "Organ"
rpm_log <- log10(rpm + 1)
cca_result <- cca(rpm_log ~ ., data = metadata)
if (is.null(cca_result$CCA) || ncol(cca_result$CCA$u) < 1) {
  stop("CCA did not produce a constrained axis.")
}

sample_order <- order(cca_result$CCA$u[, 1])
pathogen_order <- order(cca_result$CCA$v[, 1])
heatmap_matrix <- t(rpm_log)[pathogen_order, sample_order, drop = FALSE]
metadata_ordered <- metadata[sample_order, , drop = FALSE]
taxonomy_ordered <- taxonomy[pathogen_order, , drop = FALSE]

ordered_organs <- metadata_ordered$Organ
organ_gaps <- which(head(ordered_organs, -1) != tail(ordered_organs, -1))
expected_organ_gaps <- as.integer(c(42, 76, 110, 134, 154, 170, 191))
if (!identical(as.integer(organ_gaps), expected_organ_gaps)) {
  warning("Updated organ gap positions: ", paste(organ_gaps, collapse = ", "))
}

# These separators reproduce the original CCA-ranked visual blocks.
row_gaps <- c(16, 31, 48, 63)
row_gaps <- row_gaps[row_gaps < nrow(heatmap_matrix)]

taxonomy_colors <- c(
  "Acinetobacter" = "#00441B",
  "Aliivibrio" = "#238B45",
  "Amnoonviridae" = "#99D8C9",
  "Apicomplexa" = "#4D004B",
  "Nanhypoviridae" = "#88419D",
  "Arthropoda" = "#8C96C6",
  "Astroviridae" = "#BFD3E6",
  "Bornaviridae" = "#08306B",
  "Caliciviridae" = "#2171B5",
  "Circoviridae" = "#6BAED6",
  "Clostridium" = "#C6DBEF",
  "Filoviridae" = "#67001F",
  "Hepaciviridae" = "#CE1256",
  "Hepadnaviridae" = "#DF65B0",
  "Iridoviridae" = "#D4B9DA",
  "Microsporidia" = "#CC4C02",
  "Nematoda" = "#EC7014",
  "Nodaviridae" = "#FE9929",
  "Papillomaviridae" = "#FEC44F",
  "Paramyxoviridae" = "#B15928",
  "Parvoviridae" = "#8F7700",
  "Photobacterium" = "#374E55",
  "Picornaviridae" = "#8DD3C7",
  "Platyhelminthes" = "#BEBADA",
  "Pseudomonas" = "#FF69B4",
  "Retroviridae" = "#32CD32",
  "Rhabdoviridae" = "#FB8072",
  "Spinareoviridae" = "#FCCDE5",
  "Vibrio" = "#00CED1"
)

organ_colors <- c(
  "Intestine" = "#C51B7D",
  "Skin" = "#DE77AE",
  "Spleen" = "#C6DBEF",
  "Liver" = "#4292C6",
  "Muscles" = "#7FBC41",
  "Brain" = "#4D9221",
  "Gill" = "#C2A5CF",
  "Kidney" = "#762A83"
)

missing_taxonomy_colors <- setdiff(unique(taxonomy_ordered$Taxonomy), names(taxonomy_colors))
missing_organ_colors <- setdiff(unique(metadata_ordered$Organ), names(organ_colors))
if (length(missing_taxonomy_colors)) {
  stop("Missing taxonomy colors: ", paste(missing_taxonomy_colors, collapse = ", "))
}
if (length(missing_organ_colors)) {
  stop("Missing organ colors: ", paste(missing_organ_colors, collapse = ", "))
}

annotation_colors <- list(
  Taxonomy = taxonomy_colors,
  Organ = organ_colors
)

heatmap_plot <- pheatmap(
  heatmap_matrix,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = colorRampPalette(c("#F2F2F2", "#DC0000"))(100),
  annotation_col = metadata_ordered,
  annotation_row = taxonomy_ordered,
  annotation_colors = annotation_colors,
  annotation_names_row = FALSE,
  gaps_col = organ_gaps,
  gaps_row = row_gaps,
  border_color = NA,
  show_colnames = FALSE,
  fontsize = 8,
  fontsize_row = 8,
  fontfamily = "sans",
  silent = TRUE
)

row_name_index <- which(heatmap_plot$gtable$layout$name == "row_names")
if (length(row_name_index) == 1) {
  heatmap_plot$gtable$grobs[[row_name_index]]$gp$fontface <- "italic"
}

draw_heatmap <- function() {
  grid::grid.newpage()
  grid::grid.draw(heatmap_plot$gtable)
  grid::grid.text(
    "a",
    x = grid::unit(0.012, "npc"),
    y = grid::unit(0.985, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontfamily = "sans", fontsize = 18, fontface = "bold")
  )
}

grDevices::cairo_pdf(output_pdf_device, width = 10, height = 8, family = "sans")
draw_heatmap()
grDevices::dev.off()

grDevices::png(
  output_png_device,
  width = 10,
  height = 8,
  units = "in",
  res = 600,
  type = "cairo",
  bg = "white"
)
draw_heatmap()
grDevices::dev.off()
