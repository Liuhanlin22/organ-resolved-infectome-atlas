#!/usr/bin/env Rscript

# Figure 4E — pathogen richness across eight paired organs
#
# Design:
#   - Retain only individuals represented by all eight organs.
#   - Richness = number of pathogen taxa with RPM > 0 in each library.
#   - Overall test: Friedman test (eight repeated measurements per individual).
#   - Post-hoc tests: paired Wilcoxon signed-rank tests with BH correction.
#
# The script exports the extracted 896-library abundance matrix, plotting data,
# statistical results, and Figure 4E. Existing source files are not overwritten.

required_packages <- c("openxlsx", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    "\nInstall with: install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()

abundance_file <- file.path(script_dir, "Abundance_matrix_with_libs.xlsx")

if (!file.exists(abundance_file)) {
  stop(
    "Cannot find the abundance matrix:\n", abundance_file,
    "\nPlace Abundance_matrix_with_libs.xlsx beside this R script before running."
  )
}

organ_map <- c(
  B = "Brain",
  G = "Gill",
  H = "Liver",
  I = "Intestine",
  K = "Kidney",
  M = "Muscles",
  P = "Skin",
  S = "Spleen"
)

# Required display order, retained from the manuscript figure.
tissue_order <- c(
  "Intestine", "Spleen", "Skin", "Kidney",
  "Liver", "Gill", "Muscles", "Brain"
)

tissue_fill <- c(
  Intestine = "#ABD6AB",
  Spleen    = "#F2DCDB",
  Skin      = "#8ED0F4",
  Kidney    = "#EAE4AD",
  Liver     = "#F2C3A0",
  Gill      = "#CCC1DA",
  Muscles   = "#B26A62",
  Brain     = "#E2B7B9"
)

# ----------------------------- Data extraction -----------------------------

raw_data <- read.xlsx(abundance_file, sheet = 1, check.names = FALSE)
if (ncol(raw_data) < 2 || nrow(raw_data) < 1) {
  stop("The abundance workbook does not contain a valid pathogen-by-library matrix.")
}

pathogen_col <- names(raw_data)[1]
pathogen_names <- trimws(as.character(raw_data[[1]]))
if (anyNA(pathogen_names) || any(pathogen_names == "") || anyDuplicated(pathogen_names)) {
  stop("The first column must contain unique, non-empty pathogen names.")
}

library_ids <- names(raw_data)[-1]
organ_code <- sub("^.*-([BGHIKMPS])$", "\\1", library_ids)
valid_suffix <- grepl("-[BGHIKMPS]$", library_ids)
if (!all(valid_suffix)) {
  stop(
    "These library IDs do not end in -B/-G/-H/-I/-K/-M/-P/-S: ",
    paste(library_ids[!valid_suffix], collapse = ", ")
  )
}

individual_id <- sub("-[BGHIKMPS]$", "", library_ids)
library_meta <- data.frame(
  lib_id = library_ids,
  individual_id = individual_id,
  organ_code = organ_code,
  Tissue = unname(organ_map[organ_code]),
  stringsAsFactors = FALSE
)

# One and only one library per organ is required for a complete individual.
individual_qc <- aggregate(
  organ_code ~ individual_id,
  data = library_meta,
  FUN = function(x) length(unique(x))
)
names(individual_qc)[2] <- "n_distinct_organs"

individual_lib_count <- aggregate(
  lib_id ~ individual_id,
  data = library_meta,
  FUN = length
)
names(individual_lib_count)[2] <- "n_libraries"
individual_qc <- merge(individual_qc, individual_lib_count, by = "individual_id")

complete_individuals <- individual_qc$individual_id[
  individual_qc$n_distinct_organs == 8 & individual_qc$n_libraries == 8
]
complete_meta <- library_meta[library_meta$individual_id %in% complete_individuals, ]

# Stable order: individual, followed by B/G/H/I/K/M/P/S.
complete_meta$organ_code <- factor(complete_meta$organ_code, levels = names(organ_map))
complete_meta <- complete_meta[order(complete_meta$individual_id, complete_meta$organ_code), ]
complete_meta$organ_code <- as.character(complete_meta$organ_code)

if (length(complete_individuals) != 112L || nrow(complete_meta) != 896L) {
  stop(
    "Expected 112 complete individuals and 896 libraries, but found ",
    length(complete_individuals), " individuals and ", nrow(complete_meta), " libraries."
  )
}

rpm <- as.matrix(raw_data[, -1, drop = FALSE])
suppressWarnings(storage.mode(rpm) <- "numeric")
if (anyNA(rpm)) {
  stop("The abundance matrix contains missing or non-numeric RPM values.")
}
if (any(rpm < 0)) {
  stop("The abundance matrix contains negative RPM values.")
}
rownames(rpm) <- pathogen_names
colnames(rpm) <- library_ids

complete_rpm <- rpm[, complete_meta$lib_id, drop = FALSE]
richness <- colSums(complete_rpm > 0)

plot_data <- complete_meta
plot_data$Richness <- as.integer(richness[plot_data$lib_id])
plot_data$Tissue <- factor(plot_data$Tissue, levels = tissue_order)

# Confirm the paired block is exactly 112 x 8.
pair_table <- table(plot_data$individual_id, plot_data$Tissue)
if (!all(dim(pair_table) == c(112L, 8L)) || !all(pair_table == 1L)) {
  stop("The extracted data are not a complete 112-individual x 8-organ paired design.")
}

# ------------------------------ Export inputs ------------------------------

complete_matrix_export <- data.frame(
  Pathogen = pathogen_names,
  complete_rpm,
  check.names = FALSE
)

write.table(
  complete_matrix_export,
  file.path(script_dir, "Figure4E_complete_896_abundance.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

otu_export <- data.frame(
  lib_id = plot_data$lib_id,
  count = plot_data$Richness,
  stringsAsFactors = FALSE
)
write.table(
  otu_export,
  file.path(script_dir, "Figure4E_OTU_complete_896.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

annotation_export <- data.frame(
  pool = plot_data$lib_id,
  type = as.character(plot_data$Tissue),
  individual_id = plot_data$individual_id,
  stringsAsFactors = FALSE
)
write.table(
  annotation_export,
  file.path(script_dir, "Figure4E_annotation_complete_896.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ---------------------------- Paired statistics ----------------------------

friedman_result <- friedman.test(Richness ~ Tissue | individual_id, data = plot_data)
friedman_export <- data.frame(
  method = unname(friedman_result$method),
  statistic = unname(friedman_result$statistic),
  df = unname(friedman_result$parameter),
  p_value = friedman_result$p.value,
  n_individuals = length(complete_individuals),
  stringsAsFactors = FALSE
)

tissue_summary <- do.call(
  rbind,
  lapply(tissue_order, function(tissue_name) {
    x <- plot_data$Richness[plot_data$Tissue == tissue_name]
    data.frame(
      Tissue = tissue_name,
      n = length(x),
      mean = mean(x),
      sd = sd(x),
      median = median(x),
      IQR = IQR(x),
      min = min(x),
      max = max(x),
      zero_libraries = sum(x == 0),
      stringsAsFactors = FALSE
    )
  })
)

wide_richness <- reshape(
  plot_data[, c("individual_id", "Tissue", "Richness")],
  idvar = "individual_id",
  timevar = "Tissue",
  direction = "wide"
)

pair_list <- combn(tissue_order, 2, simplify = FALSE)
paired_results <- do.call(
  rbind,
  lapply(pair_list, function(pair) {
    x <- wide_richness[[paste0("Richness.", pair[1])]]
    y <- wide_richness[[paste0("Richness.", pair[2])]]
    test <- suppressWarnings(wilcox.test(
      x, y,
      paired = TRUE,
      exact = FALSE,
      correct = FALSE
    ))
    data.frame(
      group1 = pair[1],
      group2 = pair[2],
      n_pairs = length(x),
      median_difference = median(x - y),
      mean_difference = mean(x - y),
      statistic_V = unname(test$statistic),
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
  })
)
paired_results$p_adjust_BH <- p.adjust(paired_results$p_value, method = "BH")
paired_results$significance_BH <- cut(
  paired_results$p_adjust_BH,
  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
  labels = c("****", "***", "**", "*", "ns")
)

richness_export <- plot_data[, c(
  "lib_id", "individual_id", "organ_code", "Tissue", "Richness"
)]
richness_export$Tissue <- as.character(richness_export$Tissue)

write.table(
  richness_export,
  file.path(script_dir, "Figure4E_richness_complete_896.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.xlsx(
  list(
    complete_896_abundance = complete_matrix_export,
    richness_896 = richness_export,
    tissue_summary = tissue_summary,
    Friedman_test = friedman_export,
    paired_Wilcoxon_BH = paired_results,
    individual_QC = individual_qc
  ),
  file.path(script_dir, "Figure4E_complete_8_organs_richness.xlsx"),
  overwrite = TRUE
)

# --------------------------------- Figure ---------------------------------

set.seed(20260814)
y_max <- max(plot_data$Richness)

# To keep the narrow panel readable, display the six strongest significant
# paired comparisons (BH-adjusted P < 0.05). The Excel output retains all 28
# pairwise tests. Increase this value if more comparison lines are required.
n_pairwise_annotations <- 6L
significant_annotations <- paired_results[
  !is.na(paired_results$p_adjust_BH) & paired_results$p_adjust_BH < 0.05,
]
significant_annotations <- significant_annotations[
  order(significant_annotations$p_adjust_BH),
]
significant_annotations <- head(significant_annotations, n_pairwise_annotations)

if (nrow(significant_annotations) > 0) {
  significant_annotations$x_start <- match(
    significant_annotations$group1, tissue_order
  )
  significant_annotations$x_end <- match(
    significant_annotations$group2, tissue_order
  )
  significant_annotations$x_mid <- (
    significant_annotations$x_start + significant_annotations$x_end
  ) / 2

  # The smallest adjusted P value is placed at the highest level, matching
  # the visual hierarchy in the reference figure.
  annotation_step <- 0.55
  annotation_base <- y_max + 0.48
  significant_annotations$annotation_y <- annotation_base +
    rev(seq(0, by = annotation_step, length.out = nrow(significant_annotations)))
  significant_annotations$p_label <- paste0(
    "p = ",
    formatC(
      significant_annotations$p_adjust_BH,
      format = "e",
      digits = 1
    )
  )
  plot_y_max <- max(significant_annotations$annotation_y) + 0.50
} else {
  plot_y_max <- y_max + 0.55
}

p <- ggplot(plot_data, aes(x = Tissue, y = Richness, fill = Tissue)) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    adjust = 1.15,
    width = 0.78,
    alpha = 0.40,
    colour = "grey35",
    linewidth = 0.35
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.42
  ) +
  geom_point(
    position = position_jitter(width = 0.13, height = 0, seed = 20260814),
    shape = 21,
    size = 3.45,
    stroke = 0.20,
    alpha = 0.55,
    colour = "grey25"
  ) +
  geom_segment(
    data = significant_annotations,
    aes(
      x = x_start,
      xend = x_end,
      y = annotation_y,
      yend = annotation_y
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.45,
    lineend = "butt"
  ) +
  geom_text(
    data = significant_annotations,
    aes(
      x = x_mid,
      y = annotation_y + 0.11,
      label = p_label
    ),
    inherit.aes = FALSE,
    family = "Arial",
    size = 3.2,
    colour = "black",
    vjust = 0
  ) +
  scale_fill_manual(values = tissue_fill, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, plot_y_max),
    breaks = seq(0, y_max, by = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Pathogen species richness"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 35, hjust = 1, vjust = 1,
      colour = "black", size = 10
    ),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.title.y = element_text(colour = "black", size = 11.5),
    axis.line = element_line(colour = "black", linewidth = 0.45),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    plot.margin = margin(t = 9, r = 8, b = 8, l = 14)
  )

figure_width_in <- 4.4
figure_height_in <- 6.2

ggsave(
  file.path(script_dir, "Figure4E.pdf"),
  p,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(script_dir, "Figure4E.png"),
  p,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  dpi = 600,
  bg = "white"
)
