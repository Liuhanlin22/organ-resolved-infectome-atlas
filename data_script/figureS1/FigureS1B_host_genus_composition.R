#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_paths <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) "" else as.character(frame$ofile)
}, character(1))
source_paths <- source_paths[nzchar(source_paths)]
script_dir <- if (length(script_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/"))
} else if (length(source_paths) > 0) {
  dirname(normalizePath(tail(source_paths, 1), winslash = "/"))
} else {
  getwd()
}

input_file <- file.path(script_dir, "FigureS1B_host_genus_composition_input.xlsx")
output_pdf <- file.path(script_dir, "FigureS1B_host_genus_composition.pdf")
output_png <- file.path(script_dir, "FigureS1B_host_genus_composition_600dpi.png")

data <- read_excel(input_file, sheet = "main")
stopifnot(names(data)[1] == "order_order_level_lineage")
if (any(!is.finite(as.matrix(data[-1]))) || any(as.matrix(data[-1]) < 0) ||
    any(abs(rowSums(as.matrix(data[-1])) - 1) > 1e-10)) {
  stop("Host-genus proportions must be nonnegative and sum to one within each order.")
}

long_data <- data %>%
  pivot_longer(
    cols = -order_order_level_lineage,
    names_to = "host_genus",
    values_to = "proportion"
  ) %>%
  filter(proportion > 0)

order_levels <- data$order_order_level_lineage
genus_levels <- names(data)[-1]
long_data <- long_data %>%
  mutate(
    order_order_level_lineage = factor(order_order_level_lineage, levels = order_levels),
    host_genus = factor(host_genus, levels = genus_levels)
  )

host_genus_colors <- c(
  "Arius" = "#FFEBEE",
  "Benthosema" = "#FFCDD2",
  "Bregmaceros" = "#EF9A9A",
  "Chelon" = "#26C6DA",
  "Clupanodon" = "#EF5350",
  "Coilia" = "#E53935",
  "Cynoglossus" = "#FCE4EC",
  "Decapterus" = "#F8BBD0",
  "Deveximentum" = "#8E24AA",
  "Eleutheronema" = "#F48FB1",
  "Encrasicholina" = "#EC407A",
  "Gymnothorax" = "#D81B60",
  "Harpadon" = "#AD1457",
  "Ilisha" = "#880E4F",
  "Inegocia" = "#C62828",
  "Johnius" = "#E1BEE7",
  "Kumococius" = "#CE93D8",
  "Lagocephalus" = "#BA68C8",
  "Larimichthys" = "#AB47BC",
  "Lepturacanthus" = "#6A1B9A",
  "Lutjanus" = "#4A148C",
  "Moolgarda" = "#1976D2",
  "Nematalosa" = "#E3F2FD",
  "Nemipterus" = "#BBDEFB",
  "Odontamblyopus" = "#2196F3",
  "Otolithes" = "#0D47A1",
  "Pampus" = "#E0F7FA",
  "Paramonacanthus" = "#B2EBF2",
  "Pennahia" = "#4DD0E1",
  "Polydactylus" = "#0097A7",
  "Psenopsis" = "#006064",
  "Sardinella" = "#E0F2F1",
  "Saurida" = "#4DB6AC",
  "Scatophagus" = "#009688",
  "Scoliodon" = "#00796B",
  "Setipinna" = "#F9FBE7",
  "Siganus" = "#E6EE9C",
  "Sillago" = "#C0CA33",
  "Telatrygon" = "#E65100",
  "Thryssa" = "#FFF8E1",
  "Zebrias" = "#FFC107"
)

missing_colors <- setdiff(genus_levels, names(host_genus_colors))
if (length(missing_colors) > 0) {
  stop("Missing colors for: ", paste(missing_colors, collapse = ", "))
}

plot <- ggplot(
  long_data,
  aes(x = order_order_level_lineage, y = proportion, fill = host_genus)
) +
  geom_col(width = 0.70, color = NA) +
  scale_fill_manual(values = host_genus_colors, breaks = genus_levels, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = sprintf("%.2f", seq(0, 1, 0.25)),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(drop = FALSE, expand = expansion(add = 0.45)) +
  labs(x = NULL, y = "Sample composition", fill = "Host genus", tag = "b") +
  guides(
    fill = guide_legend(
      nrow = 4,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0,
      keyheight = grid::unit(0.30, "cm"),
      keywidth = grid::unit(0.32, "cm")
    )
  ) +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(
    axis.title.y = element_text(size = 12, margin = margin(r = 8)),
    axis.text.y = element_text(size = 9.5, color = "black"),
    axis.text.x = element_text(
      size = 9.5, color = "black", angle = 90, hjust = 1, vjust = 0.5
    ),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(color = "#888888", fill = NA, linewidth = 0.45),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 8.2, face = "italic"),
    legend.spacing.y = grid::unit(0.03, "cm"),
    legend.box.margin = margin(0, 0, 6, 0),
    plot.tag = element_text(size = 18, face = "bold", family = "Arial"),
    plot.tag.position = "topleft",
    plot.margin = margin(8, 8, 8, 12)
  ) +
  coord_cartesian(ylim = c(0, 1.045), clip = "off")

grDevices::pdfFonts(Arial = grDevices::pdfFonts("ArialMT")[[1]])
ggsave(output_pdf, plot = plot, width = 16.0, height = 7.0, units = "in",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(output_png, plot = plot, width = 16.0, height = 7.0, units = "in", dpi = 600, bg = "white")
