#!/usr/bin/env Rscript

# Figure 5C — number of distinct pathogens detected per fish

required <- c("openxlsx", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else normalizePath(getwd(), winslash = "/")

input_file <- file.path(script_dir, "Figure5C_input.xlsx")
if (!file.exists(input_file)) stop("Input not found: ", input_file)

pie_data <- read.xlsx(input_file, sheet = "pie_summary")
pie_data$fish_count <- as.numeric(pie_data$fish_count)
pie_data$proportion <- as.numeric(pie_data$proportion)

category_order <- c(
  "pathogen 0", "pathogen 1", "pathogen 2", "pathogen 3",
  "pathogen 4", "pathogen 5", "pathogen > 5"
)
pie_data$category <- factor(pie_data$category, levels = category_order)

if (sum(pie_data$fish_count) != 112L) stop("Figure 5C counts must sum to 112 fish.")

display_names <- c(
  "pathogen 0" = "No pathogens\ndetected",
  "pathogen 1" = "Single pathogen",
  "pathogen 2" = "Two pathogens",
  "pathogen 3" = "Three pathogens",
  "pathogen 4" = "Four pathogens",
  "pathogen 5" = "Five pathogens",
  "pathogen > 5" = "> 5 pathogens"
)

pie_colors <- c(
  "pathogen 0" = "#efedf5",
  "pathogen 1" = "#dadaeb",
  "pathogen 2" = "#bcbddc",
  "pathogen 3" = "#9e9ac8",
  "pathogen 4" = "#807dba",
  "pathogen 5" = "#6a51a3",
  "pathogen > 5" = "#54278f"
)

pie_data$label <- paste0(
  display_names[as.character(pie_data$category)], "\n(",
  pie_data$fish_count, ", ", sprintf("%.2f%%", 100 * pie_data$proportion), ")"
)

# Arrange wedges to match the reference figure: large slices on the left,
# progressively smaller slices on the right.
draw_order <- c(
  "pathogen 1", "pathogen 2", "pathogen 3", "pathogen 5",
  "pathogen 4", "pathogen > 5", "pathogen 0"
)
pie_data$draw_order <- match(as.character(pie_data$category), draw_order)
pie_data <- pie_data[order(pie_data$draw_order), ]

start_angle <- 193 * pi / 180
pie_data$angle_width <- 2 * pi * pie_data$proportion
pie_data$angle_start <- start_angle - c(0, head(cumsum(pie_data$angle_width), -1))
pie_data$angle_end <- pie_data$angle_start - pie_data$angle_width
pie_data$angle_mid <- (pie_data$angle_start + pie_data$angle_end) / 2

wedge_points <- do.call(rbind, lapply(seq_len(nrow(pie_data)), function(i) {
  theta <- seq(pie_data$angle_start[i], pie_data$angle_end[i], length.out = 120)
  data.frame(
    category = pie_data$category[i],
    x = c(0, cos(theta), 0),
    y = c(0, sin(theta), 0)
  )
}))

inside <- pie_data[as.character(pie_data$category) %in%
                     c("pathogen 0", "pathogen 1", "pathogen 2"), ]
inside$x <- 0.59 * cos(inside$angle_mid)
inside$y <- 0.59 * sin(inside$angle_mid)

outside <- pie_data[!as.character(pie_data$category) %in%
                      c("pathogen 0", "pathogen 1", "pathogen 2"), ]
outside$x0 <- 0.92 * cos(outside$angle_mid)
outside$y0 <- 0.92 * sin(outside$angle_mid)
outside$x1 <- 1.08
outside$label_x <- 1.18
outside$label_y <- c(0.55, 0.18, -0.20, -0.56)

p <- ggplot() +
  geom_polygon(
    data = wedge_points,
    aes(x = x, y = y, group = category, fill = category),
    color = "white", linewidth = 1.1
  ) +
  geom_text(
    data = inside,
    aes(x = x, y = y, label = label),
    family = "Arial", size = 3.25, lineheight = 0.92
  ) +
  geom_segment(
    data = outside,
    aes(x = x0, y = y0, xend = x1, yend = label_y),
    linewidth = 0.4, color = "black"
  ) +
  geom_segment(
    data = outside,
    aes(x = x1, y = label_y, xend = label_x - 0.03, yend = label_y),
    linewidth = 0.4, color = "black"
  ) +
  geom_text(
    data = outside,
    aes(x = label_x, y = label_y, label = label),
    hjust = 0, family = "Arial", size = 2.75, lineheight = 0.90
  ) +
  scale_fill_manual(values = pie_colors, drop = FALSE) +
  coord_fixed(xlim = c(-1.25, 2.08), ylim = c(-1.18, 1.18), clip = "off") +
  labs(tag = "c") +
  theme_void(base_family = "Arial", base_size = 10) +
  theme(
    legend.position = "none",
    plot.tag = element_text(family = "Arial", face = "bold", size = 16),
    plot.tag.position = c(0.02, 0.98),
    plot.margin = margin(8, 15, 8, 8)
  )

ggsave(file.path(script_dir, "Figure5C.pdf"), p,
       width = 6.2, height = 4.4, units = "in",
       device = grDevices::cairo_pdf, bg = "white")
ggsave(file.path(script_dir, "Figure5C.png"), p,
       width = 6.2, height = 4.4, units = "in",
       dpi = 600, bg = "white")
