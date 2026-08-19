#!/usr/bin/env Rscript

# Figure 4D — organ-specific pathogen prevalence and pathogen count
#
# Input: Abundance_matrix_with_libs.xlsx (pathogen RPM matrix)
# Only individuals represented by all eight organ libraries are retained.
#
# Organ suffixes:
#   B = Brain, G = Gill, H = Liver, I = Intestine,
#   K = Kidney, M = Muscles, P = Skin, S = Spleen
#
# Definitions:
#   Prevalence (%) = individuals with >=1 pathogen (RPM > 0) / complete
#                    individuals * 100, calculated separately for each organ.
#   Pathogen count = number of distinct pathogen rows with RPM > 0 in at least
#                    one complete individual for the indicated organ.

required_packages <- c("openxlsx", "writexl", "ggplot2")
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
    library(writexl)
    library(ggplot2)
  })
)

# ---------------------------- adjustable parameters -------------------------
font_family <- "Arial"
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
}

figure_width_in <- 4.5
figure_height_in <- 5.2
bar_width <- 0.72
line_width <- 1.15
point_size <- 3.2
line_color <- "#F8CBAD"
point_color <- "#C00000"

plot_tissue_order <- c(
  "Intestine", "Spleen", "Skin", "Kidney",
  "Liver", "Gill", "Muscles", "Brain"
)

tissue_colors <- c(
  "Intestine" = "#ABD6AB",
  "Spleen" = "#F2DCDB",
  "Skin" = "#8ED0F4",
  "Kidney" = "#EAE4AD",
  "Liver" = "#F2C3A0",
  "Gill" = "#CCC1DA",
  "Muscles" = "#B26A62",
  "Brain" = "#E2B7B9"
)

organ_code_to_tissue <- c(
  "B" = "Brain",
  "G" = "Gill",
  "H" = "Liver",
  "I" = "Intestine",
  "K" = "Kidney",
  "M" = "Muscles",
  "P" = "Skin",
  "S" = "Spleen"
)
required_organ_codes <- names(organ_code_to_tissue)

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
input_file <- file.path(script_dir, "Abundance_matrix_with_libs.xlsx")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

# ---------------------------- read and validate -----------------------------
abundance <- openxlsx::read.xlsx(
  input_file,
  sheet = 1,
  check.names = FALSE
)
if (ncol(abundance) < 2) stop("The abundance matrix has fewer than two columns.")
names(abundance)[1] <- "Pathogen"
abundance$Pathogen <- trimws(as.character(abundance$Pathogen))

if (anyDuplicated(abundance$Pathogen)) {
  stop("The abundance matrix contains duplicated pathogen names.")
}

library_ids <- names(abundance)[-1]
library_pattern <- "^(.*)-([BGHIKMPS])$"
valid_library_id <- grepl(library_pattern, library_ids)
if (any(!valid_library_id)) {
  stop(
    "Library name(s) do not end in a recognized organ suffix: ",
    paste(library_ids[!valid_library_id], collapse = ", ")
  )
}

individual_ids <- sub(library_pattern, "\\1", library_ids)
organ_codes <- sub(library_pattern, "\\2", library_ids)
library_manifest <- data.frame(
  individual_id = individual_ids,
  organ_code = organ_codes,
  tissue = unname(organ_code_to_tissue[organ_codes]),
  lib_index = library_ids,
  stringsAsFactors = FALSE
)

if (anyDuplicated(library_manifest$lib_index)) {
  stop("Duplicated library IDs were detected.")
}
if (anyDuplicated(library_manifest[, c("individual_id", "organ_code")])) {
  stop("At least one individual has duplicated libraries for the same organ.")
}

organ_table <- table(
  library_manifest$individual_id,
  factor(library_manifest$organ_code, levels = required_organ_codes)
)
organ_present <- organ_table > 0
complete_individuals <- rownames(organ_present)[rowSums(organ_present) == 8]
incomplete_individuals <- setdiff(rownames(organ_present), complete_individuals)

if (length(complete_individuals) == 0) stop("No complete eight-organ individuals were found.")

incomplete_report <- data.frame(
  individual_id = incomplete_individuals,
  organ_count = rowSums(organ_present[incomplete_individuals, , drop = FALSE]),
  missing_organ_codes = vapply(
    incomplete_individuals,
    function(id) {
      paste(required_organ_codes[!organ_present[id, ]], collapse = ",")
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)
incomplete_report$missing_tissues <- vapply(
  strsplit(incomplete_report$missing_organ_codes, ",", fixed = TRUE),
  function(x) paste(unname(organ_code_to_tissue[x]), collapse = ", "),
  character(1)
)
incomplete_report <- incomplete_report[
  order(-incomplete_report$organ_count, incomplete_report$individual_id),
  ,
  drop = FALSE
]

# Convert the pathogen x library matrix to numeric only after all names have
# been captured and validated.
rpm_matrix <- as.matrix(abundance[, -1, drop = FALSE])
suppressWarnings(storage.mode(rpm_matrix) <- "numeric")
if (anyNA(rpm_matrix)) {
  stop("The abundance matrix contains NA or non-numeric RPM values.")
}
rownames(rpm_matrix) <- abundance$Pathogen
colnames(rpm_matrix) <- library_ids

# Canonical export order: alphabetical individual ID, then B/G/H/I/K/M/P/S.
complete_individuals <- sort(complete_individuals)
complete_library_ids <- unlist(
  lapply(
    complete_individuals,
    function(id) paste0(id, "-", required_organ_codes)
  ),
  use.names = FALSE
)
if (any(!complete_library_ids %in% colnames(rpm_matrix))) {
  stop("Internal error: at least one expected complete library is missing.")
}
complete_rpm_matrix <- rpm_matrix[, complete_library_ids, drop = FALSE]

complete_manifest <- library_manifest[
  match(complete_library_ids, library_manifest$lib_index),
  ,
  drop = FALSE
]

# ---------------------------- summary statistics ----------------------------
summary_data <- do.call(
  rbind,
  lapply(plot_tissue_order, function(tissue_name) {
    organ_code <- names(organ_code_to_tissue)[
      organ_code_to_tissue == tissue_name
    ]
    tissue_library_ids <- paste0(complete_individuals, "-", organ_code)
    tissue_matrix <- rpm_matrix[, tissue_library_ids, drop = FALSE]

    positive_individual <- colSums(tissue_matrix > 0) > 0
    pathogen_detected <- rowSums(tissue_matrix > 0) > 0

    data.frame(
      Tissue = tissue_name,
      Complete_individuals = length(complete_individuals),
      Positive_individuals = sum(positive_individual),
      Prevalence_percent = sum(positive_individual) /
        length(complete_individuals) * 100,
      Pathogen_count = sum(pathogen_detected),
      stringsAsFactors = FALSE
    )
  })
)

summary_data$Tissue <- factor(summary_data$Tissue, levels = plot_tissue_order)
summary_data$Tissue_index <- seq_len(nrow(summary_data))

# Count individuals with any detected pathogen across all eight organs. All-zero
# individuals remain in the prevalence denominator, as required for an unbiased
# complete-case comparison.
individual_any_positive <- vapply(
  complete_individuals,
  function(id) {
    any(rpm_matrix[, paste0(id, "-", required_organ_codes), drop = FALSE] > 0)
  },
  logical(1)
)
individual_qc <- data.frame(
  individual_id = complete_individuals,
  Any_pathogen_detected = individual_any_positive,
  stringsAsFactors = FALSE
)

# ---------------------------- data exports ----------------------------------
summary_export <- summary_data
summary_export$Tissue <- as.character(summary_export$Tissue)
summary_export$Tissue_index <- NULL

complete_matrix_export <- data.frame(
  Pathogen = rownames(complete_rpm_matrix),
  complete_rpm_matrix,
  check.names = FALSE
)

writexl::write_xlsx(
  list(
    organ_summary = summary_export,
    complete_individuals = individual_qc,
    complete_library_manifest = complete_manifest,
    incomplete_individuals = incomplete_report,
    complete_8_organs_RPM = complete_matrix_export
  ),
  file.path(script_dir, "Figure4D_complete_8_organs_analysis.xlsx")
)

write.table(
  summary_export,
  file.path(script_dir, "Figure4D_prevalence_pathogen_count.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
  fileEncoding = "UTF-8"
)
write.table(
  complete_manifest,
  file.path(script_dir, "Figure4D_complete_library_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
  fileEncoding = "UTF-8"
)
write.table(
  incomplete_report,
  file.path(script_dir, "Figure4D_incomplete_individuals.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
  fileEncoding = "UTF-8"
)

# ---------------------------- visualization ---------------------------------
# Natural cubic interpolation is used only to draw a visually gentle line;
# red points always show the exact integer pathogen counts.
smooth_x <- seq(1, nrow(summary_data), length.out = 400)
smooth_y <- stats::splinefun(
  summary_data$Tissue_index,
  summary_data$Pathogen_count,
  method = "natural"
)(smooth_x)
smooth_line <- data.frame(
  Tissue_index = smooth_x,
  Pathogen_count = pmax(0, smooth_y)
)

y_upper <- ceiling(max(
  summary_data$Prevalence_percent,
  summary_data$Pathogen_count
) / 5) * 5
y_upper <- max(40, y_upper)

p <- ggplot(summary_data, aes(x = Tissue_index)) +
  geom_col(
    aes(y = Prevalence_percent, fill = Tissue),
    width = bar_width,
    color = NA
  ) +
  geom_line(
    data = smooth_line,
    aes(x = Tissue_index, y = Pathogen_count),
    inherit.aes = FALSE,
    color = line_color,
    linewidth = line_width,
    lineend = "round"
  ) +
  geom_point(
    aes(y = Pathogen_count),
    color = point_color,
    fill = point_color,
    shape = 21,
    size = point_size,
    stroke = 0.25
  ) +
  scale_fill_manual(values = tissue_colors, guide = "none") +
  scale_x_continuous(
    breaks = seq_along(plot_tissue_order),
    labels = plot_tissue_order,
    limits = c(0.45, length(plot_tissue_order) + 0.55),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Pathogen prevalence (%)",
    limits = c(0, y_upper),
    breaks = seq(0, y_upper, by = 5),
    expand = c(0, 0),
    sec.axis = sec_axis(~ ., name = "Pathogen count")
  ) +
  annotate(
    "rect",
    xmin = 5.35, xmax = 5.72, ymin = y_upper - 3.1, ymax = y_upper - 2.2,
    fill = tissue_colors[["Intestine"]], color = NA
  ) +
  annotate(
    "text",
    x = 5.84, y = y_upper - 2.65,
    label = "Prevalence",
    hjust = 0, vjust = 0.5,
    family = font_family, size = 3.2, color = "black"
  ) +
  annotate(
    "segment",
    x = 5.35, xend = 5.72,
    y = y_upper - 5.3, yend = y_upper - 5.3,
    color = line_color,
    linewidth = line_width,
    lineend = "round"
  ) +
  annotate(
    "point",
    x = 5.535, y = y_upper - 5.3,
    color = point_color,
    fill = point_color,
    shape = 21,
    size = 2.4,
    stroke = 0.2
  ) +
  annotate(
    "text",
    x = 5.84, y = y_upper - 5.3,
    label = "Pathogen count",
    hjust = 0, vjust = 0.5,
    family = font_family, size = 3.2, color = "black"
  ) +
  labs(x = NULL, tag = "d") +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = font_family, base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 38, hjust = 1, vjust = 1,
      color = "black", size = 9.5
    ),
    axis.text.y = element_text(color = "black", size = 9.5),
    axis.title.y.left = element_text(
      color = "black", size = 11,
      margin = margin(r = 7)
    ),
    axis.title.y.right = element_text(
      color = "black", size = 11,
      margin = margin(l = 7)
    ),
    axis.line = element_line(color = "#666666", linewidth = 0.45),
    axis.ticks = element_line(color = "#666666", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "#D9D9D9", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.tag = element_text(
      family = font_family, face = "bold", size = 18,
      margin = margin(r = 4)
    ),
    plot.tag.position = c(0.005, 0.995),
    plot.margin = margin(t = 9, r = 14, b = 8, l = 12),
    legend.position = "none"
  )

pdf_file <- file.path(script_dir, "Figure4D_prevalence_pathogen_count.pdf")
grDevices::cairo_pdf(
  pdf_file,
  width = figure_width_in,
  height = figure_height_in,
  family = font_family,
  bg = "white"
)
print(p)
grDevices::dev.off()

png_file <- file.path(script_dir, "Figure4D_prevalence_pathogen_count.png")
grDevices::png(
  png_file,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  res = 600,
  bg = "white"
)
print(p)
grDevices::dev.off()


