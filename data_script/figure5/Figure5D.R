#!/usr/bin/env Rscript

# Figure 5D — overlap of pathogen classes among pathogen-positive fish

required <- c("openxlsx", "VennDiagram", "grid")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(openxlsx)
  library(VennDiagram)
  library(grid)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else normalizePath(getwd(), winslash = "/")

input_file <- file.path(script_dir, "Figure5D_input.xlsx")
if (!file.exists(input_file)) stop("Input not found: ", input_file)

status <- read.xlsx(input_file, sheet = "class_status", check.names = FALSE)
names(status)[names(status) == "RNA.virus"] <- "RNA virus"
names(status)[names(status) == "DNA.virus"] <- "DNA virus"
class_names <- c("RNA virus", "DNA virus", "Bacteria", "Eukaryota")
status[class_names] <- lapply(status[class_names], as.numeric)

positive <- rowSums(status[class_names]) > 0
if (nrow(status) != 112L || sum(positive) != 72L) {
  stop("Expected 112 total fish and 72 pathogen-positive fish; found ",
       nrow(status), " and ", sum(positive), ".")
}

sets <- list(
  "RNA virus" = status$individual_id[status$`RNA virus` == 1],
  "DNA virus" = status$individual_id[status$`DNA virus` == 1],
  "Bacteria" = status$individual_id[status$Bacteria == 1],
  "Eukaryota" = status$individual_id[status$Eukaryota == 1]
)

venn_colors <- c(
  "RNA virus" = "#eaf3e9",
  "Bacteria" = "#b7bccd",
  "Eukaryota" = "#c4dde1",
  "DNA virus" = "#929789"
)

if (requireNamespace("futile.logger", quietly = TRUE)) {
  invisible(futile.logger::flog.threshold(futile.logger::ERROR))
}

venn_grob <- venn.diagram(
  x = sets,
  filename = NULL,
  category.names = names(sets),
  fill = unname(venn_colors[names(sets)]),
  alpha = 0.72,
  col = "#666666",
  lwd = 0.6,
  cex = 0.92,
  fontfamily = "Arial",
  cat.cex = 1.15,
  cat.fontfamily = "Arial",
  cat.fontface = "plain",
  margin = 0.07,
  print.mode = c("raw", "percent"),
  sigdigs = 2,
  disable.logging = TRUE
)

draw_figure <- function() {
  grid.newpage()
  grid.draw(venn_grob)
  grid.text("d", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
            just = c("left", "top"),
            gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 16))
}

grDevices::cairo_pdf(file.path(script_dir, "Figure5D.pdf"), width = 6.0, height = 4.7)
draw_figure()
dev.off()

png(file.path(script_dir, "Figure5D.png"), width = 3600, height = 2820,
    res = 600, type = "cairo", bg = "white")
draw_figure()
dev.off()
