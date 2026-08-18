# Figure S4C: pathogen abundance across eight paired tissue types

required_packages <- c("ggplot2", "readxl", "dplyr", "tidyr", "writexl", "ragg")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)

script_arguments <- commandArgs(trailingOnly = FALSE)
script_file_argument <- grep("^--file=", script_arguments, value = TRUE)
output_dir <- if (length(script_file_argument) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_argument[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}
abundance_file <- file.path(output_dir, "Abundance_matrix_with_libs_complete8_896.xlsx")
metadata_file <- file.path(output_dir, "Sample_metadata_complete8_896.xlsx")

tissue_order <- c(
  "Intestine", "Spleen", "Skin", "Kidney",
  "Liver", "Gill", "Muscles", "Brain"
)

tissue_colors <- c(
  Intestine = "#ABD6AB",
  Spleen = "#F2DCDB",
  Skin = "#8ED0F4",
  Kidney = "#EAE4AD",
  Liver = "#F2C3A0",
  Gill = "#CCC1DA",
  Muscles = "#B26A62",
  Brain = "#E2B7B9"
)

if (!file.exists(abundance_file) || !file.exists(metadata_file)) {
  stop("Both input workbooks must be stored beside this R script.")
}

# The abundance workbook is a 75-pathogen x 896-library RPM matrix.
abundance_matrix <- read_excel(abundance_file, sheet = 1)
pathogen_column <- names(abundance_matrix)[1]
sample_columns <- trimws(names(abundance_matrix)[-1])
names(abundance_matrix) <- c(pathogen_column, sample_columns)

if (nrow(abundance_matrix) != 75 || length(sample_columns) != 896) {
  stop("The abundance matrix must contain 75 pathogen rows and 896 sample columns.")
}
if (anyDuplicated(sample_columns)) {
  stop("Duplicated sample IDs were found in the abundance matrix.")
}

abundance_numeric <- abundance_matrix %>%
  select(all_of(sample_columns)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

if (anyNA(abundance_numeric)) {
  stop("Non-numeric or missing abundance values were found in the abundance matrix.")
}
if (any(as.matrix(abundance_numeric) < 0)) {
  stop("Negative abundance values were found in the abundance matrix.")
}

sample_abundance <- tibble(
  lib_index = sample_columns,
  abundance_RPM = colSums(as.data.frame(abundance_numeric), na.rm = FALSE)
)

metadata <- read_excel(metadata_file, sheet = "main") %>%
  transmute(
    lib_index = trimws(as.character(lib_index)),
    individual_id = trimws(as.character(individual_id)),
    tissue_types = trimws(as.character(tissue_types))
  )

if (nrow(metadata) != 896 || n_distinct(metadata$individual_id) != 112) {
  stop("The metadata must contain 896 libraries from 112 complete fish.")
}
if (anyNA(metadata) || any(metadata == "")) {
  stop("Missing sample IDs, fish IDs, or tissue labels were found in the metadata.")
}
if (anyDuplicated(metadata$lib_index)) {
  stop("Duplicated lib_index values were found in the metadata.")
}
if (!setequal(metadata$tissue_types, tissue_order)) {
  stop("The metadata tissue labels do not match the required eight tissue types.")
}

complete_check <- metadata %>%
  count(individual_id, tissue_types, name = "n_library") %>%
  complete(individual_id, tissue_types = tissue_order, fill = list(n_library = 0)) %>%
  group_by(individual_id) %>%
  summarise(
    n_library_total = sum(n_library),
    n_tissue = sum(n_library == 1),
    all_tissues_once = all(n_library == 1),
    .groups = "drop"
  )

if (any(!complete_check$all_tissues_once)) {
  stop("At least one fish does not have exactly one library for each of the eight tissues.")
}

matrix_only <- setdiff(sample_abundance$lib_index, metadata$lib_index)
metadata_only <- setdiff(metadata$lib_index, sample_abundance$lib_index)
if (length(matrix_only) > 0 || length(metadata_only) > 0) {
  stop(
    "Sample IDs do not match. Matrix-only: ", length(matrix_only),
    "; metadata-only: ", length(metadata_only)
  )
}

plot_data <- sample_abundance %>%
  inner_join(metadata, by = "lib_index") %>%
  mutate(
    tissue_types = factor(tissue_types, levels = tissue_order),
    log10_RPM_plus_1 = log10(abundance_RPM + 1)
  ) %>%
  arrange(tissue_types, individual_id)

if (nrow(plot_data) != 896 || any(table(plot_data$tissue_types) != 112)) {
  stop("Joined data must contain 112 libraries for each of the eight tissues.")
}

paired_log <- plot_data %>%
  select(individual_id, tissue_types, log10_RPM_plus_1) %>%
  pivot_wider(names_from = tissue_types, values_from = log10_RPM_plus_1) %>%
  arrange(individual_id)

paired_raw <- plot_data %>%
  select(individual_id, tissue_types, abundance_RPM) %>%
  pivot_wider(names_from = tissue_types, values_from = abundance_RPM) %>%
  arrange(individual_id)

if (nrow(paired_log) != 112 || anyNA(paired_log) || anyNA(paired_raw)) {
  stop("The paired wide-format data are incomplete.")
}

# Overall repeated-measures test across the eight tissues.
friedman_fit <- friedman.test(as.matrix(paired_log[, tissue_order]))
friedman_stats <- tibble(
  test = "Friedman rank-sum test",
  outcome = "log10(total pathogen RPM + 1)",
  n_fish = nrow(paired_log),
  n_tissues = length(tissue_order),
  chi_square = unname(friedman_fit$statistic),
  df = unname(friedman_fit$parameter),
  p_value = friedman_fit$p.value,
  kendalls_W = unname(friedman_fit$statistic) /
    (nrow(paired_log) * (length(tissue_order) - 1))
)

rank_biserial_paired <- function(difference) {
  difference <- difference[!is.na(difference) & difference != 0]
  if (length(difference) == 0) {
    return(0)
  }
  signed_ranks <- rank(abs(difference), ties.method = "average")
  w_positive <- sum(signed_ranks[difference > 0])
  w_negative <- sum(signed_ranks[difference < 0])
  (w_positive - w_negative) / (w_positive + w_negative)
}

# Two-sided paired Wilcoxon signed-rank tests; positive effects mean group1 > group2.
pair_list <- combn(tissue_order, 2, simplify = FALSE)
pairwise_stats <- bind_rows(lapply(pair_list, function(pair) {
  x <- paired_log[[pair[1]]]
  y <- paired_log[[pair[2]]]
  x_raw <- paired_raw[[pair[1]]]
  y_raw <- paired_raw[[pair[2]]]
  difference <- x - y
  fit <- wilcox.test(
    x,
    y,
    paired = TRUE,
    alternative = "two.sided",
    exact = FALSE,
    correct = TRUE
  )
  tibble(
    group1 = pair[1],
    group2 = pair[2],
    n_pairs = length(difference),
    n_zero_differences = sum(difference == 0),
    wilcoxon_V = unname(fit$statistic),
    median_difference_log10 = median(difference),
    median_difference_RPM = median(x_raw - y_raw),
    paired_rank_biserial_r = rank_biserial_paired(difference),
    p_value = fit$p.value
  )
})) %>%
  mutate(
    p_BH = p.adjust(p_value, method = "BH"),
    significant_BH = p_BH < 0.05,
    significance_BH = case_when(
      p_BH < 0.0001 ~ "****",
      p_BH < 0.001 ~ "***",
      p_BH < 0.01 ~ "**",
      p_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    effect_direction = case_when(
      paired_rank_biserial_r > 0 ~ "group1 > group2",
      paired_rank_biserial_r < 0 ~ "group1 < group2",
      TRUE ~ "no direction"
    ),
    x1 = match(group1, tissue_order),
    x2 = match(group2, tissue_order)
  )

format_p <- function(p) {
  if (p < 0.0001) {
    return(formatC(p, format = "e", digits = 1))
  }
  formatC(p, format = "f", digits = 4)
}

annotation_data <- pairwise_stats %>%
  filter(significant_BH) %>%
  mutate(
    span = x2 - x1,
    label = paste0("BH p = ", vapply(p_BH, format_p, character(1)))
  ) %>%
  arrange(span, x1, x2)

data_max <- max(plot_data$log10_RPM_plus_1)
bracket_step <- max(0.20, data_max * 0.055)
bracket_base <- data_max + 0.18
if (nrow(annotation_data) > 0) {
  annotation_data <- annotation_data %>%
    mutate(y = bracket_base + (row_number() - 1) * bracket_step)
}
y_limit <- if (nrow(annotation_data) > 0) {
  max(annotation_data$y) + 0.18
} else {
  data_max + 0.3
}

tissue_summary <- plot_data %>%
  group_by(tissue_types) %>%
  summarise(
    n_libraries = n(),
    n_fish = n_distinct(individual_id),
    zero_n = sum(abundance_RPM == 0),
    zero_percent = 100 * mean(abundance_RPM == 0),
    min_RPM = min(abundance_RPM),
    q1_RPM = quantile(abundance_RPM, 0.25),
    median_RPM = median(abundance_RPM),
    mean_RPM = mean(abundance_RPM),
    q3_RPM = quantile(abundance_RPM, 0.75),
    max_RPM = max(abundance_RPM),
    median_log10_RPM_plus_1 = median(log10_RPM_plus_1),
    .groups = "drop"
  )

figure_s4c <- ggplot(
  plot_data,
  aes(x = tissue_types, y = log10_RPM_plus_1, fill = tissue_types)
) +
  geom_violin(
    width = 0.78,
    trim = TRUE,
    scale = "width",
    alpha = 0.72,
    color = NA
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    alpha = 0.9,
    linewidth = 0.25
  ) +
  geom_point(
    position = position_jitter(width = 0.13, height = 0, seed = 20260817),
    size = 3.00,
    alpha = 0.58,
    color = "black",
    shape = 16
  ) +
  scale_fill_manual(values = tissue_colors, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, ceiling(y_limit), by = 1),
    limits = c(0, y_limit),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = expression("Pathogen abundance " * log[10] * "(RPM+1)")
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#D9D9D9", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.55),
    axis.title.y = element_text(size = 10.5, color = "black", margin = margin(r = 7)),
    axis.text.y = element_text(size = 9.5, color = "black"),
    axis.text.x = element_text(size = 11, color = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks = element_line(color = "#333333", linewidth = 0.4),
    legend.position = "none",
    plot.margin = margin(8, 8, 8, 8)
  )

if (nrow(annotation_data) > 0) {
  cap_height <- 0.045
  figure_s4c <- figure_s4c +
    geom_segment(
      data = annotation_data,
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.32,
      color = "black"
    ) +
    geom_segment(
      data = annotation_data,
      aes(x = x1, xend = x1, y = y, yend = y - cap_height),
      inherit.aes = FALSE,
      linewidth = 0.32,
      color = "black"
    ) +
    geom_segment(
      data = annotation_data,
      aes(x = x2, xend = x2, y = y, yend = y - cap_height),
      inherit.aes = FALSE,
      linewidth = 0.32,
      color = "black"
    ) +
    geom_text(
      data = annotation_data,
      aes(x = (x1 + x2) / 2, y = y + 0.035, label = label),
      inherit.aes = FALSE,
      size = 2.55,
      vjust = 0,
      color = "black"
    )
}

pdf_file <- file.path(output_dir, "FigureS4C_pathogen_abundance_by_tissue.pdf")
svg_file <- file.path(output_dir, "FigureS4C_pathogen_abundance_by_tissue.svg")
png_file <- file.path(output_dir, "FigureS4C_pathogen_abundance_by_tissue.png")
tiff_file <- file.path(output_dir, "FigureS4C_pathogen_abundance_by_tissue.tiff")

ggsave(pdf_file, figure_s4c, width = 4.35, height = 6.2, units = "in", device = cairo_pdf)
svg(svg_file, width = 4.35, height = 6.2, bg = "white", onefile = TRUE)
print(figure_s4c)
dev.off()
# Some Windows raster devices cannot write directly to a non-ASCII path. Use
# the ASCII temporary directory for rasterization, then copy the finished files
# back to the requested figure directory.
raster_dir <- file.path(tempdir(), "figureS4C_raster")
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
png_stage <- file.path(raster_dir, basename(png_file))
tiff_stage <- file.path(raster_dir, basename(tiff_file))
ggsave(
  png_stage,
  figure_s4c,
  width = 4.35,
  height = 6.2,
  units = "in",
  dpi = 600,
  bg = "white",
  device = ragg::agg_png
)
ggsave(
  tiff_stage,
  figure_s4c,
  width = 4.35,
  height = 6.2,
  units = "in",
  dpi = 600,
  bg = "white",
  compression = "lzw",
  device = "tiff"
)
if (!file.copy(png_stage, png_file, overwrite = TRUE)) {
  stop("Could not copy the PNG output to the figure directory.")
}
if (!file.copy(tiff_stage, tiff_file, overwrite = TRUE)) {
  stop("Could not copy the TIFF output to the figure directory.")
}

write.csv(
  plot_data,
  file.path(output_dir, "FigureS4C_plot_data.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  tissue_summary,
  file.path(output_dir, "FigureS4C_tissue_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  friedman_stats,
  file.path(output_dir, "FigureS4C_friedman_test.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  pairwise_stats,
  file.path(output_dir, "FigureS4C_pairwise_paired_wilcoxon_BH.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

workbook_file <- file.path(output_dir, "FigureS4C_data_and_statistics.xlsx")
writexl::write_xlsx(
  list(
    plot_data = as.data.frame(plot_data),
    tissue_summary = as.data.frame(tissue_summary),
    friedman_test = as.data.frame(friedman_stats),
    paired_wilcoxon_BH = as.data.frame(pairwise_stats),
    plotted_annotations = as.data.frame(annotation_data),
    cohort_check = as.data.frame(complete_check)
  ),
  path = workbook_file,
  col_names = TRUE,
  format_headers = TRUE
)
