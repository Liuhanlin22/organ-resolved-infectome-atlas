#!/usr/bin/env Rscript
required <- c("openxlsx", "grid")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(openxlsx); library(grid)})

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(arg) == 1) dirname(normalizePath(sub("^--file=", "", arg), winslash = "/")) else normalizePath(getwd(), winslash = "/")
input_file <- file.path(script_dir, "FigureS6A_input_updated.xlsx")
out_base <- file.path(script_dir, "FigureS6A_heatmap")

matrix_df <- read.xlsx(input_file, sheet = "abundance_RPM", check.names = FALSE)
host_df <- read.xlsx(input_file, sheet = "host_annotation", check.names = FALSE)
path_df <- read.xlsx(input_file, sheet = "pathogen_annotation", check.names = FALSE)
order_df <- read.xlsx(input_file, sheet = "host_order_colors", check.names = FALSE)

matrix_species <- trimws(gsub("\\.", " ", names(matrix_df)[-1]))
host_species <- trimws(as.character(host_df$species))
stopifnot(identical(matrix_species, host_species))
names(matrix_df)[-1] <- matrix_species
stopifnot(identical(as.character(matrix_df[[1]]), as.character(path_df$pathogen)))
stopifnot(!anyDuplicated(names(matrix_df)[-1]), !anyDuplicated(matrix_df[[1]]))
abundance <- as.matrix(matrix_df[, -1, drop = FALSE])
storage.mode(abundance) <- "double"
if (any(!is.finite(abundance)) || any(abundance < 0)) stop("abundance_RPM contains invalid values")
rownames(abundance) <- as.character(matrix_df[[1]])
colnames(abundance) <- names(matrix_df)[-1]

# The displayed scale is deliberately fixed at 0-3. Values above the upper
# display limit are clipped so the legend and the colour mapping agree.
heatmap_data <- pmin(log10(abundance + 1), 3)
colour_values <- colorRampPalette(c("#f7f7f7", "#1d4d4f"))(301)
family_palette <- c(
  Aquareoviridae = "#f4cccc", Nodaviridae = "#e69138", Acinetobacter = "#f6b26b",
  Aliivibrio = "#f9cb9c", Clostridium = "#f768a1", Photobacterium = "#93c47d",
  Vibrio = "#b2e3f1", Nematoda = "#b4a7d6", Platyhelminthes = "#fce5cd"
)
order_palette <- setNames(as.character(order_df$colour), as.character(order_df$host_order))
if (!all(unique(as.character(host_df$host_order)) %in% names(order_palette))) stop("Missing host-order colour")
if (!all(unique(as.character(path_df$family_display)) %in% names(family_palette))) stop("Missing pathogen-family colour")

draw_legend <- function(title, labels, colours, x, y_top, item_step = 0.035) {
  grid.text(title, x = unit(x, "npc"), y = unit(y_top, "npc"), just = c("left", "top"),
            gp = gpar(fontfamily = "sans", fontsize = 7.5, fontface = "bold"))
  for (i in seq_along(labels)) {
    y <- y_top - 0.027 - (i - 1) * item_step
    grid.rect(x = unit(x + 0.005, "npc"), y = unit(y, "npc"), width = unit(0.010, "npc"),
              height = unit(0.018, "npc"), gp = gpar(fill = colours[[labels[i]]], col = "#bbbbbb", lwd = 0.3))
    grid.text(labels[i], x = unit(x + 0.020, "npc"), y = unit(y, "npc"), just = c("left", "center"),
              gp = gpar(fontfamily = "sans", fontsize = 6))
  }
}

draw_figure <- function() {
  grid.newpage()
  pushViewport(viewport(x = 0, y = 0, width = 1, height = 1, just = c("left", "bottom")))
  nrow_h <- nrow(heatmap_data); ncol_h <- ncol(heatmap_data)
  left <- 0.060; bottom <- 0.090; grid_w <- 0.515; grid_h <- 0.470
  top <- bottom + grid_h; cw <- grid_w / ncol_h; ch <- grid_h / nrow_h
  family_x <- left - 0.010

  grid.text("a", x = unit(0.008, "npc"), y = unit(0.975, "npc"), just = c("left", "top"),
            gp = gpar(fontfamily = "sans", fontsize = 15, fontface = "bold"))

  # Host annotation strip and strictly vertical (90 degree) species names.
  for (j in seq_len(ncol_h)) {
    x <- left + (j - 0.5) * cw
    grid.rect(x = unit(x, "npc"), y = unit(top + 0.013, "npc"), width = unit(cw, "npc"), height = unit(0.018, "npc"),
              gp = gpar(fill = order_palette[[as.character(host_df$host_order[j])]], col = "white", lwd = 0.3))
    grid.text(colnames(heatmap_data)[j], x = unit(x, "npc"), y = unit(top + 0.038, "npc"), rot = 90,
              just = c("left", "center"), gp = gpar(fontfamily = "sans", fontsize = 6, fontface = "italic"))
    if (j > 1 && as.character(host_df$host_order[j]) != as.character(host_df$host_order[j - 1])) {
      x_line <- left + (j - 1) * cw
      grid.lines(x = unit(c(x_line, x_line), "npc"), y = unit(c(bottom, top + 0.013), "npc"),
                 gp = gpar(col = "white", lwd = 2.2))
    }
  }

  # Heatmap body, family strip, row labels, and family separators.
  for (i in seq_len(nrow_h)) {
    y <- bottom + grid_h - (i - 0.5) * ch
    fam <- as.character(path_df$family_display[i])
    grid.rect(x = unit(family_x, "npc"), y = unit(y, "npc"), width = unit(0.009, "npc"), height = unit(ch, "npc"),
              gp = gpar(fill = family_palette[[fam]], col = "white", lwd = 0.3))
    for (j in seq_len(ncol_h)) {
      value <- heatmap_data[i, j]
      idx <- min(301L, max(1L, floor(value / 3 * 300) + 1L))
      grid.rect(x = unit(left + (j - 0.5) * cw, "npc"), y = unit(y, "npc"), width = unit(cw, "npc"), height = unit(ch, "npc"),
                gp = gpar(fill = colour_values[idx], col = "white", lwd = 0.25))
    }
    grid.text(rownames(heatmap_data)[i], x = unit(left + grid_w + 0.010, "npc"), y = unit(y, "npc"),
              just = c("left", "center"), gp = gpar(fontfamily = "sans", fontsize = 6, fontface = "italic"))
    if (i > 1 && as.character(path_df$family_display[i]) != as.character(path_df$family_display[i - 1])) {
      grid.lines(x = unit(c(left, left + grid_w), "npc"), y = unit(bottom + grid_h - (i - 1) * ch, "npc"),
                 gp = gpar(col = "white", lwd = 2.2))
    }
  }
  grid.text("Host order", x = unit(family_x - 0.004, "npc"), y = unit(top + 0.013, "npc"), rot = 90,
            just = c("center", "bottom"), gp = gpar(fontfamily = "sans", fontsize = 7, fontface = "bold"))
  grid.text("Family", x = unit(family_x - 0.004, "npc"), y = unit(bottom + grid_h / 2, "npc"), rot = 90,
            just = c("center", "center"), gp = gpar(fontfamily = "sans", fontsize = 7, fontface = "bold"))

  # Legends are intentionally outside the matrix in the upper-right area.
  order_labels <- unique(as.character(host_df$host_order))
  family_labels <- unique(as.character(path_df$family_display))
  draw_legend("Host order", order_labels, order_palette, 0.680, 0.965)
  draw_legend("Pathogen taxa", family_labels, family_palette, 0.810, 0.965)

  # Fixed abundance scale with exactly the requested labels.
  bar_x <- 0.955; bar_y <- 0.575; bar_h <- 0.335; bar_w <- 0.012
  grid.text("log10(RPM+1)", x = unit(0.935, "npc"), y = unit(0.965, "npc"), just = c("left", "top"),
            gp = gpar(fontfamily = "sans", fontsize = 7.5, fontface = "bold"))
  for (k in 0:100) {
    y <- bar_y + (k + 0.5) / 101 * bar_h
    idx <- min(301L, max(1L, floor(k / 100 * 300) + 1L))
    grid.rect(x = unit(bar_x, "npc"), y = unit(y, "npc"), width = unit(bar_w, "npc"), height = unit(bar_h / 101 + 0.0005, "npc"),
              gp = gpar(fill = colour_values[idx], col = NA))
  }
  for (v in 0:3) {
    y <- bar_y + v / 3 * bar_h
    grid.lines(x = unit(c(bar_x + bar_w, bar_x + bar_w + 0.004), "npc"), y = unit(c(y, y), "npc"), gp = gpar(col = "#111111", lwd = 0.5))
    grid.text(as.character(v), x = unit(bar_x + bar_w + 0.007, "npc"), y = unit(y, "npc"), just = c("left", "center"),
              gp = gpar(fontfamily = "sans", fontsize = 6))
  }
  popViewport()
}

pdf(paste0(out_base, ".pdf"), width = 8, height = 10 / 3, useDingbats = FALSE)
draw_figure(); dev.off()
png(paste0(out_base, ".png"), width = 2400, height = 1000, res = 300, bg = "white")
draw_figure(); dev.off()
if (capabilities("cairo")) {
  tiff(paste0(out_base, ".tiff"), width = 2400, height = 1000, res = 300, compression = "lzw", type = "cairo", bg = "white")
  draw_figure(); dev.off()
}
