#!/usr/bin/env Rscript

# Figure 2B — pathogen abundance heatmap
# Input: revised Table S4 containing 51 host species, 8 tissues and 4 metadata columns.
# The script rebuilds all four annotation/order files on every run so they cannot
# become out of sync with Table S4.

required_packages <- c("openxlsx", "pheatmap", "grid")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them with install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "), "))."
  )
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(pheatmap)
  library(grid)
})

# ------------------------------ File settings ------------------------------
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
default_workbook <- "Supplement_Table4_Normalized abundance, related to figure2-revised-S3-annotations.xlsx"
input_xlsx <- if (length(args) >= 1) args[1] else file.path(script_dir, default_workbook)
output_dir <- if (length(args) >= 2) args[2] else script_dir

if (!file.exists(input_xlsx)) {
  stop(
    "Cannot find Table S4: ", input_xlsx, "\n",
    "Place the workbook beside figure2B.R or run:\n",
    "Rscript figure2B.R <TableS4.xlsx> <output_directory>"
  )
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

annotation_col_file <- file.path(output_dir, "col_heatmap_annot.txt")
annotation_row_file <- file.path(output_dir, "row_heatmap_annot.txt")
order_col_file <- file.path(output_dir, "order_data.txt")
order_row_file <- file.path(output_dir, "order_data_row.txt")
pdf_file <- file.path(output_dir, "figure2B_heatmap.pdf")
png_file <- file.path(output_dir, "figure2B_heatmap.png")

# ---------------------------- Read revised Table S4 -------------------------
raw <- openxlsx::read.xlsx(
  input_xlsx,
  sheet = 1,
  colNames = FALSE,
  skipEmptyRows = FALSE,
  skipEmptyCols = FALSE,
  check.names = FALSE
)
raw <- as.data.frame(raw, check.names = FALSE, stringsAsFactors = FALSE)

if (nrow(raw) < 5 || ncol(raw) < 10) stop("Table S4 does not have the expected structure.")
headers <- trimws(as.character(unlist(raw[4, , drop = TRUE])))
classification_col <- match("Classification", headers)
taxonomy_col <- match("Taxonomy", headers)
novelty_col <- match("Novelty", headers)
multi_col <- match("Multi-tissue distribution", headers)
if (anyNA(c(classification_col, taxonomy_col, novelty_col, multi_col))) {
  stop("Table S4 must contain Classification, Taxonomy, Novelty and Multi-tissue distribution columns.")
}

data_cols <- 2:(classification_col - 1)
pathogen_rows <- which(
  seq_len(nrow(raw)) >= 5 &
    !is.na(raw[[1]]) &
    nzchar(trimws(as.character(raw[[1]])))
)

pathogen_names <- trimws(as.character(raw[pathogen_rows, 1]))
column_ids <- headers[data_cols]
if (anyDuplicated(pathogen_names)) stop("Duplicated pathogen names were found in Table S4.")
if (anyDuplicated(column_ids)) stop("Duplicated host/tissue column names were found in Table S4.")

to_numeric <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y[is.na(y)] <- 0
  y
}
abundance_matrix <- do.call(cbind, lapply(raw[pathogen_rows, data_cols, drop = FALSE], to_numeric))
abundance_matrix <- as.matrix(abundance_matrix)
rownames(abundance_matrix) <- pathogen_names
colnames(abundance_matrix) <- column_ids

tissue_names <- c("Brain", "Gill", "Intestine", "Liver", "Kidney", "Spleen", "Muscle", "Muscles", "Skin")
host_class <- trimws(as.character(unlist(raw[2, data_cols, drop = TRUE])))
host_order <- trimws(as.character(unlist(raw[3, data_cols, drop = TRUE])))
is_tissue <- column_ids %in% tissue_names
# Keep host taxonomy and organ annotations independent. The workbook currently
# uses "Muscles" as the column ID, but the figure legend uses "Muscle".
organ_type <- rep(NA_character_, length(column_ids))
organ_type[is_tissue] <- column_ids[is_tissue]
organ_type[organ_type == "Muscles"] <- "Muscle"
host_class[is_tissue] <- NA_character_
host_order[is_tissue] <- NA_character_
if (sum(!is_tissue) != 51 || sum(is_tissue) != 8 || length(pathogen_names) != 75) {
  stop(
    "Unexpected Table S4 dimensions: ", length(pathogen_names), " pathogens, ",
    sum(!is_tissue), " host columns and ", sum(is_tissue), " tissue columns."
  )
}

classification <- trimws(as.character(raw[pathogen_rows, classification_col]))
taxa <- trimws(as.character(raw[pathogen_rows, taxonomy_col]))

# --------------------- Rebuild the four annotation files --------------------
column_annotation_out <- data.frame(
  lib_id = column_ids,
  `Host class` = host_class,
  `Host order` = host_order,
  `Organ type` = organ_type,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
row_annotation_out <- data.frame(
  Pathogen = pathogen_names,
  Classification = classification,
  Taxa = taxa,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.table(column_annotation_out, annotation_col_file, sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
write.table(row_annotation_out, annotation_row_file, sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
writeLines(paste(c("lib_id", column_ids), collapse = "\t"), order_col_file, useBytes = TRUE)
writeLines(pathogen_names, order_row_file, useBytes = TRUE)

# Read the files back and validate them before plotting.
desired_columns <- strsplit(readLines(order_col_file, warn = FALSE, encoding = "UTF-8")[1], "\t", fixed = TRUE)[[1]][-1]
desired_rows <- readLines(order_row_file, warn = FALSE, encoding = "UTF-8")
if (!setequal(desired_columns, colnames(abundance_matrix))) stop("order_data.txt does not match Table S4 columns.")
if (!setequal(desired_rows, rownames(abundance_matrix))) stop("order_data_row.txt does not match Table S4 pathogens.")

abundance_matrix <- abundance_matrix[desired_rows, desired_columns, drop = FALSE]
heatmap_data <- log10(abundance_matrix + 1)

annotation_col <- read.delim(
  annotation_col_file, header = TRUE, sep = "\t", row.names = 1,
  check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)
annotation_row <- read.delim(
  annotation_row_file, header = TRUE, sep = "\t", row.names = 1,
  check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)
annotation_col <- annotation_col[colnames(heatmap_data), c("Host class", "Host order", "Organ type"), drop = FALSE]
annotation_row <- annotation_row[rownames(heatmap_data), c("Classification", "Taxa"), drop = FALSE]

# --------------------------- Annotation palettes ----------------------------
# Same 20-order palette used in Figure 1; Clupeiformes is #1C69AF.
order_colors <- c(
  "Anguilliformes" = "#54278F",
  "Aulopiformes" = "#756BB1",
  "Carangaria incertae sedis" = "#9E9AC8",
  "Carangiformes" = "#BCBDDC",
  "Chaetodontiformes" = "#DADAEB",
  "Clupeiformes" = "#1C69AF",
  "Eupercaria incertae sedis" = "#3182BD",
  "Gadiformes" = "#6BAED6",
  "Gobiiformes" = "#9ECAE1",
  "Lutjaniformes" = "#C6DBEF",
  "Mugiliformes" = "#006D2C",
  "Myctophiformes" = "#31A354",
  "Perciformes" = "#74C476",
  "Pleuronectiformes" = "#A1D99B",
  "Scombriformes" = "#C7E9C0",
  "Siluriformes" = "#C60505",
  "Spariformes" = "#DE2D26",
  "Tetraodontiformes" = "#FC9272",
  "Carcharhiniformes" = "#FCBBA1",
  "Myliobatiformes" = "#F4D4C9"
)
class_colors <- c(
  "Actinopteri" = "#9E9AC8",
  "Chondrichthyes" = "#D8A19A"
)
organ_colors <- c(
  "Brain" = "#F57C6E",
  "Gill" = "#F2B56F",
  "Intestine" = "#FAE69E",
  "Liver" = "#84C3B7",
  "Kidney" = "#88D8DB",
  "Spleen" = "#71B7ED",
  "Muscle" = "#B8AEEB",
  "Skin" = "#F2A7DA"
)
classification_colors <- c(
  "RNA virus" = "#FBB4AE",
  "DNA virus" = "#FFD9A8",
  "Bacteria" = "#B3CDE4",
  "Eukaryota" = "#CCEAC4"
)
taxa_colors <- c(
  "Astro-Poty" = "#F44336",
  "Bunya-Arena" = "#E06666",
  "Flavi" = "#EA9999",
  "Mono-Chu" = "#F4CCCC",
  "Nido" = "#CC4C02",
  "Orthomyxo" = "#E69138",
  "Picorna" = "#F6B26B",
  "Reo" = "#F9CB9C",
  "Tombus-Noda" = "#FCE5CD",
  "Ortervirales" = "#F768A1",
  "Blubervirales" = "#8FCE00",
  "Cirlivirales" = "#93C47D",
  "Herpesvirales" = "#B6D7A8",
  "Piccovirales" = "#D9EAD3",
  "Pimascovirales" = "#225EA8",
  "Zurhausenvirales" = "#2986CC",
  "Bacillota" = "#B2E3F1",
  "Pseudomonadota" = "#7FCDBB",
  "Apicomplexa" = "#6A329F",
  "Microsporidia" = "#B4A7D6",
  "Nematoda" = "#A64D79",
  "Platyhelminthes" = "#C27BA0",
  "Arthropoda" = "#EAD1DC"
)

non_missing_levels <- function(x) unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
missing_order_colors <- setdiff(non_missing_levels(annotation_col[["Host order"]]), names(order_colors))
missing_class_colors <- setdiff(non_missing_levels(annotation_col[["Host class"]]), names(class_colors))
missing_organ_colors <- setdiff(non_missing_levels(annotation_col[["Organ type"]]), names(organ_colors))
missing_taxa_colors <- setdiff(unique(annotation_row$Taxa), names(taxa_colors))
missing_classification_colors <- setdiff(unique(annotation_row$Classification), names(classification_colors))
if (length(missing_order_colors) > 0) stop("Missing Host order colors: ", paste(missing_order_colors, collapse = ", "))
if (length(missing_class_colors) > 0) stop("Missing Host class colors: ", paste(missing_class_colors, collapse = ", "))
if (length(missing_organ_colors) > 0) stop("Missing Organ type colors: ", paste(missing_organ_colors, collapse = ", "))
if (length(missing_taxa_colors) > 0) stop("Missing Taxa colors: ", paste(missing_taxa_colors, collapse = ", "))
if (length(missing_classification_colors) > 0) stop("Missing Classification colors: ", paste(missing_classification_colors, collapse = ", "))

annotation_col[["Host class"]] <- factor(annotation_col[["Host class"]], levels = c("Actinopteri", "Chondrichthyes"))
annotation_col[["Host order"]] <- factor(annotation_col[["Host order"]], levels = unique(host_order[!is.na(host_order)]))
annotation_col[["Organ type"]] <- factor(
  annotation_col[["Organ type"]],
  levels = c("Brain", "Gill", "Intestine", "Liver", "Kidney", "Spleen", "Muscle", "Skin")
)
annotation_row$Classification <- factor(annotation_row$Classification, levels = c("RNA virus", "DNA virus", "Bacteria", "Eukaryota"))
annotation_row$Taxa <- factor(annotation_row$Taxa, levels = unique(taxa))

annotation_colors <- list()
annotation_colors[["Host class"]] <- class_colors
annotation_colors[["Host order"]] <- order_colors
annotation_colors[["Organ type"]] <- organ_colors
annotation_colors[["Classification"]] <- classification_colors
annotation_colors[["Taxa"]] <- taxa_colors

# Major row gaps: between four Classification groups and between eukaryotic Taxa.
classification_chr <- as.character(annotation_row$Classification)
taxa_chr <- as.character(annotation_row$Taxa)
row_transition <- which(classification_chr[-length(classification_chr)] != classification_chr[-1])
eukaryote_transition <- which(
  classification_chr[-length(classification_chr)] == "Eukaryota" &
    classification_chr[-1] == "Eukaryota" &
    taxa_chr[-length(taxa_chr)] != taxa_chr[-1]
)
gaps_row <- sort(unique(c(row_transition, eukaryote_transition)))

# Column gaps: between host orders, and between individual tissue columns.
order_chr <- ifelse(is_tissue, "Tissue", as.character(annotation_col[["Host order"]]))
order_transition <- which(order_chr[-length(order_chr)] != order_chr[-1])
tissue_positions <- which(order_chr == "Tissue")
tissue_gaps <- if (length(tissue_positions) > 1) tissue_positions[-length(tissue_positions)] else integer(0)
gaps_col <- sort(unique(c(order_transition, tissue_gaps)))

# ------------------------------ Draw heatmap -------------------------------
max_value <- max(heatmap_data, na.rm = TRUE)
if (!is.finite(max_value) || max_value <= 0) max_value <- 1
heatmap_breaks <- seq(0, max_value, length.out = 201)
legend_breaks <- pretty(c(0, max_value), n = 5)
legend_breaks <- legend_breaks[legend_breaks >= 0 & legend_breaks <= max_value]

heatmap_result <- pheatmap::pheatmap(
  heatmap_data,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("#F7F6F5", "#DD1332"))(200),
  breaks = heatmap_breaks,
  legend_breaks = legend_breaks,
  legend_labels = format(legend_breaks, trim = TRUE, scientific = FALSE),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  gaps_row = gaps_row,
  gaps_col = gaps_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  labels_row = gsub("_+", " ", rownames(heatmap_data)),
  annotation_names_col = TRUE,
  annotation_names_row = TRUE,
  annotation_legend = TRUE,
  drop_levels = FALSE,
  # Do not draw borders around individual heatmap cells. Only the group gaps
  # defined by gaps_row and gaps_col remain, matching the original figure.
  border_color = NA,
  na_col = "#FFFFFF",
  fontsize_row = 20,
  legend = TRUE,
  silent = TRUE
)

draw_heatmap <- function() {
  grid::grid.newpage()
  grid::grid.draw(heatmap_result$gtable)
  grid::grid.text(
    "b", x = grid::unit(0.006, "npc"), y = grid::unit(0.995, "npc"),
    just = c("left", "top"), gp = grid::gpar(fontfamily = "Arial", fontsize = 12, fontface = "bold")
  )
}

if (capabilities("cairo")) {
  grDevices::cairo_pdf(pdf_file, width = 22, height = 21, family = "Arial", bg = "white")
} else {
  grDevices::pdf(pdf_file, width = 22, height = 21, family = "sans", useDingbats = FALSE)
}
draw_heatmap()
grDevices::dev.off()

png_type <- if (capabilities("cairo")) "cairo" else getOption("bitmapType")
grDevices::png(png_file, width = 22, height = 21, units = "in", res = 300, bg = "white", type = png_type)
draw_heatmap()
grDevices::dev.off()
