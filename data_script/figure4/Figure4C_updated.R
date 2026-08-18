#!/usr/bin/env Rscript

# Figure 4C — brain-library pathogen heatmap (without host photographs)
#
# This script rebuilds Figure 4C directly from:
#   1. the revised Table S1 host/sample annotation;
#   2. Abundance_matrix_with_libs.xlsx.
#
# It automatically:
#   - keeps Brain libraries only;
#   - retains the 13 Figure 4C pathogens in the specified row order;
#   - keeps libraries with RPM > 0 for at least one selected pathogen;
#   - updates host species/order/class from the revised Table S1;
#   - writes the updated matrix, sorting files and annotation files;
#   - draws only the heatmap and legends (no fish photographs).

required_packages <- c("openxlsx", "writexl", "pheatmap", "gtable")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    "\nInstall them with: install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "), "))"
  )
}

suppressWarnings(
  suppressPackageStartupMessages({
    library(openxlsx)
    library(pheatmap)
  })
)

# ---------------------------- adjustable parameters -------------------------
output_font_family <- "Arial"
# Leave grid text family empty so pheatmap can calculate its layout on the
# temporary device; the Cairo PDF device below supplies Arial as the default.
font_family <- ""
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
}

figure_width_in <- 12.0
figure_height_in <- 7.0
row_font_size <- 9.0
column_font_size <- 8.5
general_font_size <- 9.0
column_label_angle <- "45"

heatmap_colors <- grDevices::colorRampPalette(
  c("#F7F7F7", "#D8D0E5", "#9E8AC5", "#542788")
)(200)

# Colors are consistent with the revised host-order palette used elsewhere.
host_order_colors <- c(
  "Anguilliformes" = "#54278F",
  "Aulopiformes" = "#756BB1",
  "Clupeiformes" = "#1C69AF",
  "Eupercaria incertae sedis" = "#3182BD",
  "Gadiformes" = "#6BAED6",
  "Gobiiformes" = "#9ECAE1",
  "Myctophiformes" = "#31A354",
  "Perciformes" = "#74C476",
  "Siluriformes" = "#C60505",
  "Tetraodontiformes" = "#FC9272"
)

host_class_colors <- c(
  "Actinopteri" = "#B9A8C2",
  "Chondrichthyes" = "#B0BBCB"
)

pathogen_taxa_colors <- c(
  "Bornaviridae" = "#F4CCCC",
  "Filoviridae" = "#EA9999",
  "Amnoonviridae" = "#E69138",
  "Caliciviridae" = "#F6B26B",
  "Picornaviridae" = "#FFD9A8",
  "Spinareoviridae" = "#F9CB9C",
  "Nodaviridae" = "#FCE5CD",
  "Retroviridae" = "#F768A1",
  "Nenyaviridae" = "#93C47D",
  "Clostridiaceae" = "#B2E3F1",
  "Microsporidia" = "#B4A7D6"
)

classification_colors <- c(
  "RNA virus" = "#FBB4AE",
  "DNA virus" = "#FFD9A8",
  "Bacteria" = "#B3CDE4",
  "Eukaryota" = "#CCEAC4"
)

# This is the original Figure 4C row order, updated to the current abundance
# matrix nomenclature (Pseudoloma sp. RWXLT -> Pseudoloma sp. FJ).
selected_pathogens <- c(
  "Bombay duck fish bornavirus",
  "Nanhai filovirus 2",
  "Nanhai amnoonvirus 1",
  "Arius arius calicivirus",
  "Gymnothorax reevesii picornavirus",
  "Nanhai reo-like virus 4",
  "Nanhai reo-like virus 7",
  "Asian seabass Nervous Necrosis Virus",
  "Nanhai retrovirus 1",
  "Nanhai retrovirus 2",
  "Nanhai circo-like virus 2",
  "Clostridium sp. Y",
  "Pseudoloma sp. FJ"
)

# These annotations were checked against the revised Table S3 taxonomy.
pathogen_taxa_map <- c(
  "Bombay duck fish bornavirus" = "Bornaviridae",
  "Nanhai filovirus 2" = "Filoviridae",
  "Nanhai amnoonvirus 1" = "Amnoonviridae",
  "Arius arius calicivirus" = "Caliciviridae",
  "Gymnothorax reevesii picornavirus" = "Picornaviridae",
  "Nanhai reo-like virus 4" = "Spinareoviridae",
  "Nanhai reo-like virus 7" = "Spinareoviridae",
  "Asian seabass Nervous Necrosis Virus" = "Nodaviridae",
  "Nanhai retrovirus 1" = "Retroviridae",
  "Nanhai retrovirus 2" = "Retroviridae",
  "Nanhai circo-like virus 2" = "Nenyaviridae",
  "Clostridium sp. Y" = "Clostridiaceae",
  "Pseudoloma sp. FJ" = "Microsporidia"
)

classification_map <- c(
  "Bombay duck fish bornavirus" = "RNA virus",
  "Nanhai filovirus 2" = "RNA virus",
  "Nanhai amnoonvirus 1" = "RNA virus",
  "Arius arius calicivirus" = "RNA virus",
  "Gymnothorax reevesii picornavirus" = "RNA virus",
  "Nanhai reo-like virus 4" = "RNA virus",
  "Nanhai reo-like virus 7" = "RNA virus",
  "Asian seabass Nervous Necrosis Virus" = "RNA virus",
  "Nanhai retrovirus 1" = "RNA virus",
  "Nanhai retrovirus 2" = "RNA virus",
  "Nanhai circo-like virus 2" = "DNA virus",
  "Clostridium sp. Y" = "Bacteria",
  "Pseudoloma sp. FJ" = "Eukaryota"
)

# ---------------------------- path handling ---------------------------------
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    if (!identical(script_path, "-") && file.exists(script_path)) {
      return(dirname(normalizePath(script_path, winslash = "/")))
    }
  }
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
table_s1_file <- file.path(
  script_dir,
  "Supplement_Table1_Info_sample, related to figure1-revised.xlsx"
)
abundance_file <- file.path(script_dir, "Abundance_matrix_with_libs.xlsx")
if (!file.exists(table_s1_file) || !file.exists(abundance_file)) {
  stop("Required input files must be stored beside this R script.")
}

# ---------------------------- input and validation --------------------------
# Table S1 has one title row; its actual header is Excel row 2.
s1_raw <- openxlsx::read.xlsx(
  table_s1_file,
  sheet = 1,
  startRow = 2,
  check.names = FALSE
)
if (ncol(s1_raw) < 11) stop("The revised Table S1 has fewer than 11 columns.")

# Use stable column positions so punctuation normalization by openxlsx cannot
# break the mapping.
sample_info <- data.frame(
  lib_index = trimws(as.character(s1_raw[[2]])),
  tissue_types = trimws(as.character(s1_raw[[3]])),
  Host_class = trimws(as.character(s1_raw[[6]])),
  Host_order = trimws(as.character(s1_raw[[7]])),
  family = trimws(as.character(s1_raw[[8]])),
  genus = trimws(as.character(s1_raw[[9]])),
  species = trimws(as.character(s1_raw[[10]])),
  species_abbreviation = trimws(as.character(s1_raw[[11]])),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (anyNA(sample_info$lib_index) || any(!nzchar(sample_info$lib_index))) {
  stop("Table S1 contains an empty lib_index.")
}
if (anyDuplicated(sample_info$lib_index)) {
  stop("Table S1 contains duplicated lib_index values.")
}

abundance <- openxlsx::read.xlsx(
  abundance_file,
  sheet = 1,
  check.names = FALSE
)
if (ncol(abundance) < 2) stop("The abundance matrix has fewer than two columns.")
names(abundance)[1] <- "Pathogen"
abundance$Pathogen <- trimws(as.character(abundance$Pathogen))

if (anyDuplicated(abundance$Pathogen)) {
  stop("The abundance matrix contains duplicated pathogen names.")
}
missing_pathogens <- setdiff(selected_pathogens, abundance$Pathogen)
if (length(missing_pathogens) > 0) {
  stop(
    "Selected pathogen(s) missing from abundance matrix: ",
    paste(missing_pathogens, collapse = ", ")
  )
}

brain_info <- sample_info[
  tolower(sample_info$tissue_types) == "brain" &
    sample_info$lib_index %in% names(abundance),
  ,
  drop = FALSE
]
if (nrow(brain_info) == 0) stop("No Brain libraries could be matched to the abundance matrix.")

row_index <- match(selected_pathogens, abundance$Pathogen)
brain_col_index <- match(brain_info$lib_index, names(abundance))
brain_matrix <- as.matrix(
  abundance[row_index, brain_col_index, drop = FALSE]
)
suppressWarnings(storage.mode(brain_matrix) <- "numeric")
if (anyNA(brain_matrix)) {
  stop("Selected abundance values contain NA or non-numeric entries.")
}
rownames(brain_matrix) <- selected_pathogens
colnames(brain_matrix) <- brain_info$lib_index

# Reproduce the original Figure 4C selection criterion: retain a Brain library
# if at least one of the 13 selected pathogens has RPM > 0.
keep_library <- colSums(brain_matrix > 0) > 0
rpm_matrix <- brain_matrix[, keep_library, drop = FALSE]
selected_info <- brain_info[match(colnames(rpm_matrix), brain_info$lib_index), , drop = FALSE]

if (ncol(rpm_matrix) == 0) stop("No positive Brain libraries remain after filtering.")

# Class -> order -> species -> individual -> lib_index ordering. This keeps
# Actinopteri first and Chondrichthyes last if both occur in future updates.
class_order <- c("Actinopteri", "Chondrichthyes")
class_rank <- match(selected_info$Host_class, class_order)
class_rank[is.na(class_rank)] <- length(class_order) + 1L
individual_number <- suppressWarnings(as.integer(
  sub(".*-([0-9]+)-B$", "\\1", selected_info$lib_index)
))
individual_number[is.na(individual_number)] <- 999999L

sort_index <- order(
  class_rank,
  selected_info$Host_order,
  selected_info$species,
  individual_number,
  selected_info$lib_index
)
selected_info <- selected_info[sort_index, , drop = FALSE]
rpm_matrix <- rpm_matrix[, selected_info$lib_index, drop = FALSE]
individual_number <- individual_number[sort_index]

# Append the individual number only when more than one selected library belongs
# to the same host species, matching the old Figure 4C naming convention.
species_frequency <- table(selected_info$species)
selected_info$Display_label <- ifelse(
  species_frequency[selected_info$species] > 1,
  paste(selected_info$species, individual_number),
  selected_info$species
)
if (anyDuplicated(selected_info$Display_label)) {
  stop("Display labels are not unique after adding individual numbers.")
}
colnames(rpm_matrix) <- selected_info$Display_label

# ---------------------------- annotations -----------------------------------
annotation_col <- data.frame(
  `Host order` = selected_info$Host_order,
  `Host class` = selected_info$Host_class,
  row.names = selected_info$Display_label,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

annotation_row <- data.frame(
  `Pathogen taxa` = unname(pathogen_taxa_map[selected_pathogens]),
  Classification = unname(classification_map[selected_pathogens]),
  row.names = selected_pathogens,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (anyNA(annotation_col) || anyNA(annotation_row)) {
  stop("Generated annotations contain NA values.")
}

used_orders <- unique(annotation_col[["Host order"]])
missing_order_colors <- setdiff(used_orders, names(host_order_colors))
if (length(missing_order_colors) > 0) {
  stop(
    "Missing color(s) for revised host order(s): ",
    paste(missing_order_colors, collapse = ", ")
  )
}

used_classes <- unique(annotation_col[["Host class"]])
missing_class_colors <- setdiff(used_classes, names(host_class_colors))
if (length(missing_class_colors) > 0) {
  stop(
    "Missing color(s) for revised host class(es): ",
    paste(missing_class_colors, collapse = ", ")
  )
}

annotation_colors <- list(
  `Host order` = host_order_colors[used_orders],
  `Host class` = host_class_colors[used_classes],
  `Pathogen taxa` = pathogen_taxa_colors[
    unique(annotation_row[["Pathogen taxa"]])
  ],
  Classification = classification_colors[
    unique(annotation_row[["Classification"]])
  ]
)

heatmap_matrix <- log10(rpm_matrix + 1)

column_group_lengths <- rle(annotation_col[["Host order"]])$lengths
gaps_col <- if (length(column_group_lengths) > 1) {
  cumsum(column_group_lengths)[-length(column_group_lengths)]
} else {
  NULL
}

row_group_lengths <- rle(annotation_row[["Pathogen taxa"]])$lengths
gaps_row <- if (length(row_group_lengths) > 1) {
  cumsum(row_group_lengths)[-length(row_group_lengths)]
} else {
  NULL
}

# ---------------------------- exports ---------------------------------------
matrix_export <- data.frame(
  Pathogen = rownames(rpm_matrix),
  rpm_matrix,
  check.names = FALSE
)
log_matrix_export <- data.frame(
  Pathogen = rownames(heatmap_matrix),
  heatmap_matrix,
  check.names = FALSE
)
sample_export <- selected_info[, c(
  "lib_index", "Display_label", "tissue_types", "Host_class", "Host_order",
  "family", "genus", "species", "species_abbreviation"
)]
pathogen_export <- data.frame(
  Pathogen = rownames(annotation_row),
  annotation_row,
  check.names = FALSE,
  row.names = NULL
)

writexl::write_xlsx(
  list(
    selected_RPM = matrix_export,
    selected_log10_RPM_plus1 = log_matrix_export,
    selected_sample_metadata = sample_export,
    selected_pathogen_annotation = pathogen_export
  ),
  file.path(script_dir, "Figure4C_viral_data_pheatmap_updated.xlsx")
)

writeLines(
  paste(c("Pathogen", selected_info$Display_label), collapse = "\t"),
  file.path(script_dir, "Figure4C_order_data_updated.txt"),
  useBytes = TRUE
)
writeLines(
  selected_pathogens,
  file.path(script_dir, "Figure4C_order_data_row_updated.txt"),
  useBytes = TRUE
)
write.table(
  data.frame(sample_label = rownames(annotation_col), annotation_col, check.names = FALSE),
  file.path(script_dir, "Figure4C_col_heatmap_annot_updated.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
  fileEncoding = "UTF-8"
)
write.table(
  pathogen_export,
  file.path(script_dir, "Figure4C_row_heatmap_annot_updated.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
  fileEncoding = "UTF-8"
)

positive_pairs <- which(rpm_matrix > 0, arr.ind = TRUE)
verification <- data.frame(
  Pathogen = rownames(rpm_matrix)[positive_pairs[, "row"]],
  Display_label = colnames(rpm_matrix)[positive_pairs[, "col"]],
  lib_index = selected_info$lib_index[positive_pairs[, "col"]],
  Host_order = selected_info$Host_order[positive_pairs[, "col"]],
  RPM = rpm_matrix[positive_pairs],
  stringsAsFactors = FALSE
)
write.table(
  verification,
  file.path(script_dir, "Figure4C_positive_library_check.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
  fileEncoding = "UTF-8"
)

# ---------------------------- heatmap ---------------------------------------
make_heatmap <- function() {
  hm <- pheatmap::pheatmap(
    heatmap_matrix,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    color = heatmap_colors,
    border_color = NA,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    gaps_col = gaps_col,
    gaps_row = gaps_row,
    show_colnames = TRUE,
    show_rownames = TRUE,
    labels_col = colnames(heatmap_matrix),
    labels_row = rownames(heatmap_matrix),
    angle_col = column_label_angle,
    fontsize = general_font_size,
    fontsize_row = row_font_size,
    fontsize_col = column_font_size,
    legend = TRUE,
    annotation_legend = TRUE,
    annotation_names_col = TRUE,
    # Row-annotation meanings are already stated in the legends; hiding their
    # vertical names prevents clipping at the lower-left corner.
    annotation_names_row = FALSE,
    treeheight_col = 0,
    treeheight_row = 0,
    silent = TRUE,
    fontfamily = font_family
  )

  # Italicize biological names only; annotation and legend text remain plain.
  row_name_index <- which(hm$gtable$layout$name == "row_names")
  col_name_index <- which(hm$gtable$layout$name == "col_names")
  if (length(row_name_index) == 1) {
    hm$gtable$grobs[[row_name_index]]$gp$fontface <- "italic"
    hm$gtable$grobs[[row_name_index]]$gp$fontfamily <- font_family
  }
  if (length(col_name_index) == 1) {
    hm$gtable$grobs[[col_name_index]]$gp$fontface <- "italic"
    hm$gtable$grobs[[col_name_index]]$gp$fontfamily <- font_family
  }

  # pheatmap does not provide a built-in continuous-legend title. Add the
  # scientific scale name above its color bar.
  legend_index <- which(hm$gtable$layout$name == "legend")
  if (length(legend_index) == 1) {
    legend_column <- hm$gtable$layout$l[legend_index]
    hm$gtable$widths[legend_column] <- grid::unit(72, "pt")
    old_legend <- hm$gtable$grobs[[legend_index]]
    old_legend$vp <- grid::viewport(
      y = grid::unit(1, "npc") - grid::unit(15, "pt"),
      just = "top"
    )
    hm$gtable$grobs[[legend_index]] <- grid::grobTree(
      old_legend,
      grid::textGrob(
        expression(log[10](RPM + 1)),
        x = grid::unit(0, "npc"),
        y = grid::unit(1, "npc") - grid::unit(1, "pt"),
        just = c("left", "top"),
        gp = grid::gpar(
          fontsize = general_font_size,
          fontface = "bold",
          fontfamily = font_family,
          col = "black"
        )
      )
    )
  }

  # The first 45-degree column label otherwise touches the device boundary.
  hm$gtable <- gtable::gtable_add_cols(
    hm$gtable,
    grid::unit(62, "pt"),
    pos = 0
  )
  hm
}

pdf_file <- file.path(script_dir, "Figure4C_heatmap_updated.pdf")
grDevices::cairo_pdf(
  pdf_file,
  width = figure_width_in,
  height = figure_height_in,
  family = output_font_family,
  bg = "white"
)
heatmap_object <- make_heatmap()
grid::grid.newpage()
grid::grid.draw(heatmap_object$gtable)
grDevices::dev.off()

png_file <- file.path(script_dir, "Figure4C_heatmap_updated.png")
grDevices::png(
  png_file,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  res = 600,
  bg = "white"
)
heatmap_object <- make_heatmap()
grid::grid.newpage()
grid::grid.draw(heatmap_object$gtable)
grDevices::dev.off()
