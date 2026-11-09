
required_packages <- c("openxlsx", "dplyr", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
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
figure_width <- 12
figure_height <- 11.6
canvas_width <- 1328
canvas_height <- 1284
matrix_left <- 274
matrix_right <- 1166
matrix_top <- 221
matrix_bottom <- 1274
right_left <- 1178
right_right <- 1308
band_top <- 44
band_bottom <- 74
top_baseline <- 211
top_plot_top <- band_bottom
font_scale <- figure_width * 72 / canvas_width

sample_x <- matrix_left + (sample_summary$Sample_order - 0.5) *
  (matrix_right - matrix_left) / sample_count
pathogen_y <- matrix_top + (pathogen_summary$Pathogen_order - 0.5) *
  (matrix_bottom - matrix_top) / pathogen_count
row_step <- (matrix_bottom - matrix_top) / pathogen_count
column_step <- (matrix_right - matrix_left) / sample_count

figure_s5 <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    xscale = c(0, canvas_width), yscale = c(canvas_height, 0),
    gp = grid::gpar(fontfamily = "Arial"), clip = "off"
  ))
  text_at <- function(label, x, y, size, just = "centre", face = "plain", rot = 0) {
    grid::grid.text(label, x = x, y = y, default.units = "native",
      just = just, rot = rot,
      gp = grid::gpar(fontsize = size * font_scale, fontface = face, col = "black"))
  }
  rect_at <- function(left, top, width, height, fill, border = NA, line_width = 0.65) {
    grid::grid.rect(x = left, y = top, width = grid::unit(width / canvas_width, "npc"),
      height = grid::unit(height / canvas_height, "npc"),
      just = c("left", "top"), default.units = "native",
      gp = grid::gpar(fill = fill, col = border, lwd = line_width))
  }
  line_at <- function(x0, y0, x1, y1, colour = "black", width = 0.65) {
    grid::grid.segments(x0, y0, x1, y1, default.units = "native",
      gp = grid::gpar(col = colour, lwd = width))
  }
  text_at("Taxa", 46, 55, 28, "left", "bold")
  for (i in seq_along(category_order)) {
    legend_y <- 89 + (i - 1) * 33
    rect_at(47, legend_y - 10, 36, 19, category_colors[[category_order[i]]])
    text_at(category_order[i], 93, legend_y, 22, "left")
  }
  for (group in group_order) {
    positions <- sample_summary$Sample_order[sample_summary$Multiplicity_group == group]
    if (!length(positions)) next
    left <- matrix_left + (min(positions) - 1) * column_step
    right <- matrix_left + max(positions) * column_step
    center <- (left + right) / 2
    rect_at(left, band_top, right - left, band_bottom - band_top,
      group_band_colors[[group]], "black")
    if (group %in% c("Single", "Twice", "Triple")) {
      text_at(group, center, (band_top + band_bottom) / 2, 22)
    } else if (group == "Quadruple") {
      text_at(group, center - 9, 15, 16, "right")
      line_at(center - 6, 26, center, band_top + 12)
    } else {
      text_at("Quintuple and more", center + 4, 30, 16, "left")
      line_at(center + 4, 39, center - 3, band_top + 12)
    }
  }
  top_height <- top_baseline - top_plot_top
  for (i in seq_len(sample_count)) {
    height <- sample_summary$log10_RPM_plus1[i] / 4 * top_height
    rect_at(sample_x[i] - column_step * 0.37, top_baseline - height,
      column_step * 0.74, height,
      top_bar_colors[[as.character(sample_summary$Multiplicity_group[i])]])
  }
  rect_at(matrix_left, top_plot_top, matrix_right - matrix_left,
    top_height, NA, "black")
  for (tick in 0:4) {
    y <- top_baseline - tick / 4 * top_height
    line_at(matrix_left - 4, y, matrix_left, y)
    text_at(as.character(tick), matrix_left - 8, y, 18, "right")
  }
  text_at(expression(log[10](RPM + 1)), 237,
    (top_baseline + top_plot_top) / 2, 20, rot = 90)
  for (i in seq_len(pathogen_count)) {
    text_at(pathogen_summary$Display_label[i], matrix_left - 19,
      pathogen_y[i], 14.1, "right", "italic")
    rect_at(matrix_left - 14, pathogen_y[i] - row_step / 2,
      13, row_step, category_colors[[as.character(pathogen_summary$Classification[i])]])
  }
  ordered_binary <- binary[pathogen_summary$Pathogen, sample_summary$Sample, drop = FALSE]
  for (j in seq_len(sample_count)) {
    grid::grid.circle(x = rep(sample_x[j], pathogen_count), y = pathogen_y,
      r = grid::unit(3.5, "native"), default.units = "native",
      gp = grid::gpar(fill = "#D5D5D5", col = NA))
    positive <- which(ordered_binary[, j] > 0)
    if (length(positive) > 1) {
      line_at(sample_x[j], min(pathogen_y[positive]), sample_x[j],
        max(pathogen_y[positive]), "#505050", 0.8)
    }
    grid::grid.circle(x = rep(sample_x[j], length(positive)), y = pathogen_y[positive],
      r = grid::unit(5.0, "native"), default.units = "native",
      gp = grid::gpar(fill = "#333333", col = NA))
  }
  for (tick in 0:4) {
    x <- right_left + tick / 4 * (right_right - right_left)
    line_at(x, matrix_top - 2, x, matrix_bottom, "#D0D0D0", 0.5)
    line_at(x, matrix_top - 7, x, matrix_top - 2)
    text_at(as.character(tick), x, matrix_top - 17, 18)
  }
  line_at(right_left, matrix_top - 2, right_right, matrix_top - 2)
  text_at(expression(log[10](RPM + 1)), (right_left + right_right) / 2,
    matrix_top - 42, 20)
  for (i in seq_len(pathogen_count)) {
    width <- pathogen_summary$log10_RPM_plus1[i] / 4 * (right_right - right_left)
    rect_at(right_left, pathogen_y[i] - row_step * 0.36, width, row_step * 0.72,
      right_bar_colors[[as.character(pathogen_summary$Multiplicity_group[i])]])
  }
  grid::popViewport()
}

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

grDevices::cairo_pdf(paste0(output_base, ".pdf"), width = figure_width,
  height = figure_height, family = "Arial", bg = "white")
figure_s5()
grDevices::dev.off()
grDevices::png(paste0(output_base, ".png"), width = figure_width,
  height = figure_height, units = "in", res = 600, type = "cairo", bg = "white")
figure_s5()
grDevices::dev.off()
grDevices::tiff(paste0(output_base, ".tiff"), width = figure_width,
  height = figure_height, units = "in", res = 600, type = "cairo",
  compression = "lzw", bg = "white")
figure_s5()
grDevices::dev.off()
