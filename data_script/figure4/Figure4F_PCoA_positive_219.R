#!/usr/bin/env Rscript

# Figure 4F: PCoA of pathogen-positive libraries
# Repeated measures are controlled by restricting PERMANOVA permutations
# within fish identity. Output: Figure4F_PCoA_positive_219.pdf only.

required_packages <- c("openxlsx", "ggplot2", "vegan", "permute")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install missing packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
  library(vegan)
  library(permute)
})

# Use the script folder as the working folder, so the script is portable on GitHub.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}

input_file <- Sys.getenv(
  "FIGURE4F_ABUNDANCE_FILE",
  unset = file.path(script_dir, "Abundance_matrix_with_libs.xlsx")
)
output_file <- file.path(script_dir, "Figure4F_PCoA_positive_219.pdf")

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nPlace Abundance_matrix_with_libs.xlsx beside this script, ",
    "or set the FIGURE4F_ABUNDANCE_FILE environment variable."
  )
}

# ------------------------------- Settings --------------------------------

n_permutations <- 9999L
set.seed(20260814)

organ_map <- c(
  B = "Brain", G = "Gill", H = "Liver", I = "Intestine",
  K = "Kidney", M = "Muscles", P = "Skin", S = "Spleen"
)

organ_order <- c(
  "Brain", "Gill", "Intestine", "Kidney",
  "Liver", "Muscles", "Skin", "Spleen"
)

organ_colors <- c(
  Brain = "#DD5F60", Gill = "#3c6e80", Intestine = "#6a51a3",
  Kidney = "#ec7014", Liver = "#3690C0", Muscles = "#DD3497",
  Skin = "#7BCCC4", Spleen = "#C79494"
)

# -------------------------- Read and select data --------------------------

raw <- read.xlsx(input_file, sheet = 1, check.names = FALSE)
if (ncol(raw) < 2) stop("The input must contain pathogen names and library columns.")

pathogen_names <- trimws(as.character(raw[[1]]))
library_ids <- colnames(raw)[-1]

if (anyDuplicated(pathogen_names) || any(pathogen_names == "")) {
  stop("Pathogen names in the first column must be unique and non-empty.")
}
if (!all(grepl("-[BGHIKMPS]$", library_ids))) {
  stop("All library IDs must end in -B/-G/-H/-I/-K/-M/-P/-S.")
}

abundance <- as.matrix(raw[, -1, drop = FALSE])
suppressWarnings(storage.mode(abundance) <- "numeric")
if (anyNA(abundance) || any(abundance < 0)) {
  stop("The abundance matrix contains missing, non-numeric, or negative values.")
}
rownames(abundance) <- pathogen_names
colnames(abundance) <- library_ids

organ_code <- sub("^.*-([BGHIKMPS])$", "\\1", library_ids)
metadata <- data.frame(
  Sample = library_ids,
  Fish_ID = sub("-[BGHIKMPS]$", "", library_ids),
  Organ = unname(organ_map[organ_code]),
  stringsAsFactors = FALSE
)

# Retain fish represented exactly once in each of the eight organs.
fish_organ_table <- table(metadata$Fish_ID, metadata$Organ)
complete_fish <- rownames(fish_organ_table)[
  ncol(fish_organ_table) == 8 &
    rowSums(fish_organ_table == 1L) == 8L &
    rowSums(fish_organ_table) == 8L
]

metadata <- metadata[metadata$Fish_ID %in% complete_fish, , drop = FALSE]
metadata$Organ <- factor(metadata$Organ, levels = organ_order)
metadata <- metadata[order(metadata$Organ, metadata$Fish_ID), , drop = FALSE]

if (length(complete_fish) != 112L || nrow(metadata) != 896L) {
  stop(
    "Expected 112 complete fish and 896 libraries; found ",
    length(complete_fish), " fish and ", nrow(metadata), " libraries."
  )
}

library_by_pathogen <- t(abundance[, metadata$Sample, drop = FALSE])
positive <- rowSums(library_by_pathogen > 0) > 0
analysis_matrix <- library_by_pathogen[positive, , drop = FALSE]
analysis_meta <- metadata[positive, , drop = FALSE]

if (nrow(analysis_matrix) != 219L) {
  stop("Expected 219 pathogen-positive libraries; found ", nrow(analysis_matrix), ".")
}

# ---------------------------- PCoA + PERMANOVA ----------------------------

distance <- vegdist(log10(analysis_matrix + 1), method = "bray")
pcoa <- wcmdscale(distance, k = 2, eig = TRUE, add = "lingoes")

scores <- as.data.frame(pcoa$points[, 1:2, drop = FALSE])
colnames(scores) <- c("PCoA1", "PCoA2")
scores$Sample <- rownames(scores)

plot_data <- merge(analysis_meta, scores, by = "Sample", sort = FALSE)
plot_data <- plot_data[match(analysis_meta$Sample, plot_data$Sample), ]
plot_data$Organ <- factor(plot_data$Organ, levels = organ_order)

positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]
axis_percent <- 100 * pcoa$eig[1:2] / sum(positive_eigenvalues)

# Restricted permutations preserve the repeated-measures structure.
permutation_scheme <- how(nperm = n_permutations)
setBlocks(permutation_scheme) <- factor(analysis_meta$Fish_ID)

permanova <- adonis2(
  distance ~ Organ,
  data = analysis_meta,
  permutations = permutation_scheme
)

r_squared <- permanova$R2[1]
p_value <- permanova$`Pr(>F)`[1]
p_text <- if (p_value < 0.001) {
  "p < 0.001"
} else {
  paste0("p = ", format.pval(p_value, digits = 3))
}

statistics_label <- paste0(
  "PERMANOVA\n",
  "R\u00b2 = ", sprintf("%.3f", r_squared), "\n",
  p_text
)

# --------------------------------- Figure ---------------------------------

p <- ggplot(plot_data, aes(PCoA1, PCoA2, colour = Organ, fill = Organ)) +
  geom_point(shape = 21, size = 3.5, stroke = 0.30, alpha = 0.62) +
  scale_colour_manual(values = organ_colors, drop = FALSE) +
  scale_fill_manual(values = organ_colors, drop = FALSE) +
  labs(
    x = sprintf("PCoA1 (%.1f%%)", axis_percent[1]),
    y = sprintf("PCoA2 (%.1f%%)", axis_percent[2]),
    colour = "Organs",
    fill = "Organs"
  ) +
  annotate(
    "text", x = -Inf, y = Inf, label = statistics_label,
    hjust = -0.08, vjust = 1.10, family = "Arial", size = 4.2
  ) +
  coord_equal(clip = "off") +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.height = grid::unit(0.54, "cm"),
    axis.text = element_text(colour = "black", size = 9.5),
    axis.title = element_text(colour = "black", size = 11.5),
    axis.line = element_line(colour = "black", linewidth = 0.45),
    axis.ticks = element_line(colour = "black", linewidth = 0.40),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(
  output_file, p,
  width = 7.5, height = 5.6, units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

