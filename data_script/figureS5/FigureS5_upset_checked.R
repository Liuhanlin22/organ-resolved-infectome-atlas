# Figure S5: pathogen-by-individual UpSet-style abundance plot
#
# Detection is defined as RPM > 0. The input contains the 68 curated pathogen
# taxa with canonical names from Figure6A_input.xlsx and 72 pathogen-positive
# individuals from the 112 complete individuals.

required_packages <- c("openxlsx", "ggplot2", "dplyr", "tidyr", "cowplot")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
})

script_arguments <- commandArgs(trailingOnly = FALSE)
script_file_argument <- grep("^--file=", script_arguments, value = TRUE)
script_dir <- if (length(script_file_argument) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_argument[1]), winslash = "/"))
} else {
  getwd()
}

input_name <- "FigureS5_input_Figure6A_names.xlsx"
output_dir <- script_dir
input_file <- file.path(output_dir, input_name)
if (!file.exists(input_file)) {
  stop("Cannot locate ", input_name, " beside the R script.")
}

category_order <- c("RNA virus", "DNA virus", "Bacteria", "Eukaryota")
category_colors <- c(
  "RNA virus" = "#B80F4A",
  "DNA virus" = "#4D4A78",
  "Bacteria" = "#F06A21",
  "Eukaryota" = "#9DA393"
)
group_order <- c("Single", "Twice", "Triple", "Quadruple", "Quintuple and more")
top_bar_colors <- c(
  "Single" = "#FBE2D7",
  "Twice" = "#FAC0AD",
  "Triple" = "#FA987E",
  "Quadruple" = "#F87364",
  "Quintuple and more" = "#EF5B5B"
)
group_band_colors <- c(
  "Single" = "#E8E1EE",
  "Twice" = "#DCE9F0",
  "Triple" = "#E1EFD9",
  "Quadruple" = "#F7E3BE",
  "Quintuple and more" = "#F2A0B3"
)
right_bar_colors <- c(
  "Single" = "#93989F",
  "Twice" = "#AFBCCF",
  "Triple" = "#8DB5CC",
  "Quadruple" = "#67A4C3",
  "Quintuple and more" = "#35639A"
)

multiplicity_group <- function(count) {
  label <- ifelse(
    count == 1, "Single",
    ifelse(
      count == 2, "Twice",
      ifelse(count == 3, "Triple", ifelse(count == 4, "Quadruple", "Quintuple and more"))
    )
  )
  factor(label, levels = group_order)
}

abundance <- read.xlsx(input_file, sheet = "abundance", check.names = FALSE)
annotation <- read.xlsx(input_file, sheet = "pathogen_annotation", check.names = FALSE)
individual_metadata <- read.xlsx(input_file, sheet = "individual_metadata", check.names = FALSE)

required_annotation_columns <- c("Pathogen", "Display_label", "Classification")
required_individual_columns <- c(
  "Sample", "Species", "Host_label", "Host_class", "Host_order", "Family", "Genus"
)
if (names(abundance)[1] != "Pathogen") {
  stop("The abundance sheet must start with a Pathogen column.")
}
if (!all(required_annotation_columns %in% names(annotation))) {
  stop("The pathogen_annotation sheet lacks required columns.")
}
if (!all(required_individual_columns %in% names(individual_metadata))) {
  stop("The individual_metadata sheet lacks required columns.")
}

abundance$Pathogen <- trimws(as.character(abundance$Pathogen))
annotation$Pathogen <- trimws(as.character(annotation$Pathogen))
sample_columns <- names(abundance)[-1]
if (anyDuplicated(abundance$Pathogen) || anyDuplicated(sample_columns)) {
  stop("Pathogen and individual identifiers must be unique.")
}

rpm_data <- lapply(abundance[sample_columns], function(value) {
  suppressWarnings(as.numeric(value))
})
rpm <- as.matrix(as.data.frame(rpm_data, check.names = FALSE))
colnames(rpm) <- sample_columns
rownames(rpm) <- abundance$Pathogen
if (anyNA(rpm) || any(rpm < 0)) {
  stop("The abundance matrix contains missing, nonnumeric, or negative values.")
}
if (!identical(dim(rpm), c(68L, 72L))) {
  stop("Expected 68 classified pathogens and 72 pathogen-positive individuals.")
}
if (!setequal(rownames(rpm), annotation$Pathogen)) {
  stop("Pathogen annotation does not match the abundance matrix.")
}
if (!setequal(colnames(rpm), trimws(as.character(individual_metadata$Sample)))) {
  stop("Individual metadata does not match the abundance matrix.")
}
if (any(!annotation$Classification %in% category_order)) {
  stop("Unexpected pathogen classification values were found.")
}

binary <- (rpm > 0) * 1L

sample_summary <- tibble(
  Sample = colnames(rpm),
  Pathogen_count = as.integer(colSums(binary)),
  Total_RPM = as.numeric(colSums(rpm))
) %>%
  mutate(log10_RPM_plus1 = log10(Total_RPM + 1)) %>%
  arrange(Pathogen_count, desc(log10_RPM_plus1)) %>%
  mutate(
    Sample_order = row_number(),
    Multiplicity_group = multiplicity_group(Pathogen_count)
  )

pathogen_summary <- tibble(
  Pathogen = rownames(rpm),
  Host_count = as.integer(rowSums(binary)),
  Total_RPM = as.numeric(rowSums(rpm))
) %>%
  mutate(log10_RPM_plus1 = log10(Total_RPM + 1)) %>%
  left_join(
    annotation %>% select(all_of(required_annotation_columns)),
    by = "Pathogen"
  ) %>%
  arrange(Host_count, desc(log10_RPM_plus1)) %>%
  mutate(
    Pathogen_order = row_number(),
    y = n() - Pathogen_order + 1,
    Multiplicity_group = multiplicity_group(Host_count),
    Classification = factor(Classification, levels = category_order)
  )

matrix_frame <- as.data.frame(binary, check.names = FALSE)
matrix_frame$Pathogen <- rownames(binary)
matrix_long <- matrix_frame %>%
  pivot_longer(-Pathogen, names_to = "Sample", values_to = "Presence") %>%
  left_join(sample_summary %>% select(Sample, Sample_order), by = "Sample") %>%
  left_join(pathogen_summary %>% select(Pathogen, y), by = "Pathogen")

rpm_frame <- as.data.frame(rpm, check.names = FALSE)
rpm_frame$Pathogen <- rownames(rpm)
detected_pairs <- rpm_frame %>%
  pivot_longer(-Pathogen, names_to = "Sample", values_to = "RPM") %>%
  filter(RPM > 0) %>%
  left_join(
    annotation %>% select(Pathogen, Classification),
    by = "Pathogen"
  )

sample_count <- nrow(sample_summary)
pathogen_count <- nrow(pathogen_summary)
group_ranges <- sample_summary %>%
  group_by(Multiplicity_group) %>%
  summarise(
    xmin = min(Sample_order) - 0.5,
    xmax = max(Sample_order) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(
      as.character(Multiplicity_group) == "Quintuple and more",
      "Quintuple\nand more",
      as.character(Multiplicity_group)
    )
  )

group_band <- ggplot(group_ranges) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = Multiplicity_group),
    color = "black",
    linewidth = 0.35
  ) +
  geom_text(
    aes(x = (xmin + xmax) / 2, y = 0.5, label = label),
    family = "Arial",
    size = 3.0,
    lineheight = 0.85
  ) +
  scale_fill_manual(values = group_band_colors, drop = FALSE) +
  scale_x_continuous(limits = c(0.5, sample_count + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

top_plot <- ggplot(sample_summary, aes(x = Sample_order, y = log10_RPM_plus1)) +
  geom_col(aes(fill = Multiplicity_group), width = 0.76) +
  scale_fill_manual(values = top_bar_colors, drop = FALSE) +
  scale_x_continuous(limits = c(0.5, sample_count + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 4), breaks = 0:4, expand = c(0, 0)) +
  labs(x = NULL, y = expression(log[10](RPM + 1))) +
  theme_classic(base_family = "Arial", base_size = 8) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 8.5, margin = margin(r = 4)),
    axis.text.y = element_text(size = 7.5, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
    plot.margin = margin(0, 0, 0, 0)
  )

legend_data <- tibble(
  Classification = factor(category_order, levels = category_order),
  y = rev(seq_along(category_order))
)
taxa_legend <- ggplot(legend_data) +
  geom_tile(aes(x = 0.11, y = y, fill = Classification), width = 0.13, height = 0.52) +
  geom_text(
    aes(x = 0.21, y = y, label = as.character(Classification)),
    hjust = 0,
    family = "Arial",
    size = 3.3
  ) +
  annotate(
    "text",
    x = 0.02,
    y = 4.8,
    label = "Taxa",
    hjust = 0,
    vjust = 1,
    family = "Arial",
    fontface = "bold",
    size = 4.8
  ) +
  scale_fill_manual(values = category_colors, drop = FALSE) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 5), expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

label_plot <- ggplot(pathogen_summary) +
  geom_text(
    aes(x = 0.93, y = y, label = Display_label),
    hjust = 1,
    family = "Arial",
    size = 1.95
  ) +
  geom_tile(
    aes(x = 0.975, y = y, fill = Classification),
    width = 0.035,
    height = 0.84
  ) +
  scale_fill_manual(values = category_colors, drop = FALSE) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, pathogen_count + 0.5), expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

segments <- matrix_long %>%
  filter(Presence == 1) %>%
  group_by(Sample, Sample_order) %>%
  summarise(ymin = min(y), ymax = max(y), point_count = n(), .groups = "drop") %>%
  filter(point_count > 1)

matrix_plot <- ggplot(matrix_long, aes(x = Sample_order, y = y)) +
  geom_point(color = "#D5D5D5", size = 0.75) +
  geom_segment(
    data = segments,
    aes(x = Sample_order, xend = Sample_order, y = ymin, yend = ymax),
    inherit.aes = FALSE,
    color = "#4A4A4A",
    linewidth = 0.28,
    alpha = 0.82
  ) +
  geom_point(
    data = matrix_long %>% filter(Presence == 1),
    color = "#222222",
    size = 1.25
  ) +
  scale_x_continuous(limits = c(0.5, sample_count + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, pathogen_count + 0.5), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

right_plot <- ggplot(
  pathogen_summary,
  aes(x = log10_RPM_plus1, y = y, fill = Multiplicity_group)
) +
  geom_col(width = 0.72, orientation = "y") +
  scale_fill_manual(values = right_bar_colors, drop = FALSE) +
  scale_x_continuous(
    position = "top",
    limits = c(0, 4),
    breaks = 0:4,
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0.5, pathogen_count + 0.5), expand = c(0, 0)) +
  labs(x = expression(log[10](RPM + 1)), y = NULL) +
  theme_minimal(base_family = "Arial", base_size = 8) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 7.5, color = "black"),
    axis.title.x = element_text(size = 8.5, margin = margin(b = 3)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#D0D0D0", linewidth = 0.3),
    plot.margin = margin(0, 0, 0, 0)
  )

top_center <- plot_grid(
  group_band,
  top_plot,
  ncol = 1,
  rel_heights = c(0.23, 1),
  align = "v",
  axis = "lr"
)
blank_plot <- ggplot() + theme_void()
top_row <- plot_grid(
  taxa_legend,
  top_center,
  blank_plot,
  ncol = 3,
  rel_widths = c(3.35, 8, 2.05)
)
bottom_row <- plot_grid(
  label_plot,
  matrix_plot,
  right_plot,
  ncol = 3,
  rel_widths = c(3.35, 8, 2.05),
  align = "h",
  axis = "tb"
)
figure_s5 <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(1.6, 8.4)
)

write.csv(
  sample_summary,
  file.path(output_dir, "FigureS5_sample_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  pathogen_summary %>% select(-y),
  file.path(output_dir, "FigureS5_pathogen_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  detected_pairs,
  file.path(output_dir, "FigureS5_detected_pairs.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

output_base <- file.path(output_dir, "FigureS5_pathogen_individual_upset")
ggsave(
  paste0(output_base, ".pdf"),
  plot = figure_s5,
  width = 12,
  height = 10.7,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)
ggsave(
  paste0(output_base, ".png"),
  plot = figure_s5,
  width = 12,
  height = 10.7,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  paste0(output_base, ".tiff"),
  plot = figure_s5,
  width = 12,
  height = 10.7,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    paste0(output_base, ".svg"),
    plot = figure_s5,
    width = 12,
    height = 10.7,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}
