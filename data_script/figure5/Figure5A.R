#!/usr/bin/env Rscript

# Figure 5A: individual-level pathogen richness and abundance
# Input: Figure5A_input.xlsx (stored beside this script)

required_packages <- c("openxlsx", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install missing packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}

input_file <- file.path(script_dir, "Figure5A_input.xlsx")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

plot_data <- read.xlsx(input_file, sheet = "plot_data", check.names = FALSE)
brackets <- read.xlsx(input_file, sheet = "order_brackets", check.names = FALSE)

classification_order <- c("RNA virus", "DNA virus", "Bacteria", "Eukaryota")
taxa_colors <- c(
  "RNA virus" = "#EFB6C8",
  "DNA virus" = "#ADD8E6",
  "Bacteria" = "#F3C6A5",
  "Eukaryota" = "#B9DDD4"
)
abundance_color <- "#9694C4"

plot_data$classification <- factor(
  plot_data$classification,
  levels = classification_order
)
plot_data <- plot_data[order(plot_data$x_index, plot_data$classification), ]

upper_data <- plot_data
lower_data <- plot_data
individual_labels <- unique(
  plot_data[, c("x_index", "individual_id", "individual_label")]
)
individual_labels <- individual_labels[order(individual_labels$x_index), ]
individual_totals <- aggregate(pathogen_count ~ x_index, plot_data, sum)
names(individual_totals)[2] <- "bar_top"
individual_labels <- merge(
  individual_labels,
  individual_totals,
  by = "x_index",
  all.x = TRUE,
  sort = FALSE
)
individual_labels <- individual_labels[order(individual_labels$x_index), ]

max_count <- max(
  aggregate(pathogen_count ~ x_index, plot_data, sum)$pathogen_count
)
max_abundance <- max(
  aggregate(log10_group_RPM_plus1 ~ x_index, plot_data, sum)$log10_group_RPM_plus1
)

upper_axis_max <- ceiling(max_count + 1)
lower_axis_max <- ceiling(max_abundance)
bracket_y <- upper_axis_max + 4.0
label_y <- bracket_y + 0.18
plot_top <- label_y + 5.2
plot_bottom <- -lower_axis_max - 0.35

brackets$bracket_y <- bracket_y
brackets$label_y <- label_y
boundaries <- brackets[-nrow(brackets), , drop = FALSE]

p <- ggplot() +
  geom_vline(
    data = boundaries,
    aes(xintercept = boundary),
    linetype = "dashed",
    colour = "#777777",
    linewidth = 0.28
  ) +
  geom_col(
    data = upper_data,
    aes(x = x_index, y = pathogen_count, fill = classification),
    width = 0.86,
    colour = "black",
    linewidth = 0.16
  ) +
  geom_col(
    data = lower_data,
    aes(
      x = x_index,
      y = -log10_group_RPM_plus1,
      group = classification
    ),
    width = 0.86,
    fill = abundance_color,
    colour = "#7775A7",
    linewidth = 0.12
  ) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.38) +
  geom_text(
    data = individual_labels,
    aes(x = x_index, y = bar_top + 0.12, label = individual_label),
    angle = 90,
    hjust = 0,
    vjust = 0.5,
    family = "Arial",
    fontface = "italic",
    size = 2.15
  ) +
  geom_segment(
    data = brackets,
    aes(x = start - 0.36, xend = end + 0.36, y = bracket_y, yend = bracket_y),
    linewidth = 0.38,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = brackets,
    aes(
      x = start - 0.36, xend = start - 0.36,
      y = bracket_y, yend = bracket_y - 0.34
    ),
    linewidth = 0.38,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = brackets,
    aes(
      x = end + 0.36, xend = end + 0.36,
      y = bracket_y, yend = bracket_y - 0.34
    ),
    linewidth = 0.38,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = brackets,
    aes(x = mid, y = label_y, label = order),
    angle = 90,
    hjust = 0.5,
    vjust = 0.5,
    family = "Arial",
    size = 2.55,
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = -4.0,
    y = upper_axis_max / 2,
    label = "Number of pathogen species",
    angle = 90,
    family = "Arial",
    size = 3.5
  ) +
  annotate(
    "text",
    x = -4.0,
    y = -lower_axis_max / 2,
    label = "Pathogen abundance\nlog10(RPM+1)",
    angle = 90,
    family = "Arial",
    size = 3.5
  ) +
  annotate(
    "text",
    x = -3.2,
    y = plot_top,
    label = "A",
    hjust = 0,
    vjust = 1,
    family = "Arial",
    fontface = "bold",
    size = 7
  ) +
  annotate(
    "rect",
    xmin = 1.0,
    xmax = 2.4,
    ymin = -lower_axis_max + 0.28,
    ymax = -lower_axis_max + 0.68,
    fill = abundance_color,
    colour = NA
  ) +
  annotate(
    "text",
    x = 2.8,
    y = -lower_axis_max + 0.48,
    label = "Abundance of pathogens",
    hjust = 0,
    family = "Arial",
    size = 3.15
  ) +
  scale_fill_manual(
    name = "Taxa",
    values = taxa_colors,
    breaks = classification_order,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = NULL,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(plot_bottom, plot_top),
    breaks = c(
      seq(-lower_axis_max, -1, by = 1),
      0:upper_axis_max
    ),
    labels = function(x) abs(x),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    xlim = c(0.5, max(plot_data$x_index) + 0.5),
    clip = "off"
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      override.aes = list(colour = NA)
    )
  ) +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 8.5, colour = "black"),
    axis.ticks.y = element_line(linewidth = 0.35),
    axis.line.y = element_line(linewidth = 0.42),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9.5),
    legend.key.height = grid::unit(0.45, "cm"),
    legend.key.width = grid::unit(0.45, "cm"),
    plot.margin = margin(t = 34, r = 12, b = 8, l = 54)
  )

ggsave(
  file.path(script_dir, "Figure5A.pdf"),
  p,
  width = 18,
  height = 7.2,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(script_dir, "Figure5A.png"),
  p,
  width = 18,
  height = 7.2,
  units = "in",
  dpi = 300,
  bg = "white"
)
