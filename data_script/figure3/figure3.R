#!/usr/bin/env Rscript

# Figure 3 — organ-specific pathogen prevalence and abundance
#
# Inputs (in the same directory as this script):
#   figure3_brain.xlsx, figure3_gill.xlsx, figure3_intestine.xlsx,
#   figure3_kidney.xlsx, figure3_liver.xlsx, figure3_muscle.xlsx,
#   figure3_skin.xlsx, figure3_spleen.xlsx
#   Supplement_Table3_Pathogens info, related to figure2-figure3-revised.xlsx
#
# Output:
#   figure3_combined_2rows_4cols.pdf

required_packages <- c("openxlsx", "ggplot2", "gridExtra")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    "\nInstall them with: install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

# ---------------------------- user-adjustable parameters --------------------
font_family <- "Arial"

classification_colors <- c(
  "RNA virus" = "#E7B8C1",
  "DNA virus" = "#A65383",
  "Bacteria"  = "#72558C",
  "Eukaryota" = "#91B3A4"
)
abundance_color <- "#F28E3B"
known_star_color <- "#E31A1C"

species_font_size <- 7.2
axis_font_size <- 8.5
organ_title_size <- 11
bar_width <- 0.78
bar_alpha <- 1
line_width <- 0.45
point_size <- 1.35

combined_width_in <- 22.0
combined_height_in <- 15.0

# ---------------------------- path handling ---------------------------------
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
setwd(script_dir)

s3_file <- file.path(
  script_dir,
  "Supplement_Table3_Pathogens info, related to figure2-figure3-revised.xlsx"
)

organ_order <- c(
  "Brain", "Gill", "Intestine", "Kidney",
  "Liver", "Muscle", "Skin", "Spleen"
)

organ_files <- setNames(
  file.path(script_dir, paste0("figure3_", tolower(organ_order), ".xlsx")),
  organ_order
)

missing_inputs <- c(s3_file, organ_files)[!file.exists(c(s3_file, organ_files))]
if (length(missing_inputs) > 0) {
  stop("Missing input file(s):\n", paste(missing_inputs, collapse = "\n"))
}

# ---------------------------- name and data helpers -------------------------
normalize_name <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- gsub("[[:space:]]+", " ", x)
  aliases <- c(
    "Toxocara sp. CTSZ-2" = "Toxocara sp. CTSZ",
    "Polydactylus sextarius picornavirus 1" =
      "Polydactylus sextarius picornavirus",
    "Cystoisospora sp. DSLT-2" = "Cystoisospora sp. DSLT-2",
    "Sarcocystidae sp. DSLT-2" = "Cystoisospora sp. DSLT-2",
    "Photobacterium damselae QLXSD" = "Photobacterium sp. QLXSD"
  )
  hit <- match(x, names(aliases))
  x[!is.na(hit)] <- unname(aliases[hit[!is.na(hit)]])
  x
}

strip_rank_prefix <- function(x) {
  sub("^[a-z]__?", "", as.character(x))
}

normalize_header_name <- function(z) {
  z <- tolower(trimws(as.character(z)))
  z <- gsub("[^a-z0-9]+", " ", z)
  trimws(gsub("[[:space:]]+", " ", z))
}

find_summary_sheet <- function(path) {
  sheets <- getSheetNames(path)
  candidates <- character(0)
  for (sheet in sheets) {
    dat <- read.xlsx(path, sheet = sheet, rows = 1:2, colNames = TRUE)
    headers <- normalize_header_name(names(dat))
    has_count <- any(grepl("positive individuals|non zero count", headers))
    has_abundance <- any(grepl("average log10|non zero avg log", headers))
    if (has_count && has_abundance && length(headers) <= 10) {
      candidates <- c(candidates, sheet)
    }
  }
  if (length(candidates) == 0) {
    stop("No four-column summary sheet found in ", basename(path))
  }
  if ("main" %in% candidates) return("main")
  candidates[[1]]
}

find_one_column <- function(nms, pattern, path) {
  normalized <- normalize_header_name(nms)
  hit <- grep(pattern, normalized, ignore.case = TRUE)
  if (length(hit) != 1) {
    stop(
      "Expected one column matching '", pattern, "' in ", basename(path),
      "; detected columns: ", paste(nms, collapse = " | ")
    )
  }
  nms[hit]
}

nice_integer_breaks <- function(maximum, n = 4) {
  step <- max(1, ceiling(maximum / n))
  unique(c(seq(0, maximum, by = step), maximum))
}

read_s3 <- function(path) {
  # Detect the header row instead of assuming it is always row 2. This also
  # tolerates openxlsx versions that convert spaces in column names to dots.
  preview <- read.xlsx(
    path, sheet = "Sheet3", rows = 1:10, colNames = FALSE,
    skipEmptyRows = FALSE, skipEmptyCols = FALSE
  )
  header_row <- NA_integer_
  for (i in seq_len(nrow(preview))) {
    keys <- normalize_header_name(unlist(preview[i, ], use.names = FALSE))
    if ("species name" %in% keys && "genome type" %in% keys && "kingdom" %in% keys) {
      header_row <- i
      break
    }
  }
  if (is.na(header_row)) {
    stop("Cannot locate the Table S3 header row in sheet 'Sheet3'.")
  }

  x <- read.xlsx(path, sheet = "Sheet3", startRow = header_row, check.names = FALSE)
  column_keys <- normalize_header_name(names(x))
  find_s3_column <- function(key) {
    hit <- which(column_keys == key)
    if (length(hit) != 1) {
      stop(
        "Table S3 column not found or duplicated: '", key,
        "'. Detected columns: ", paste(names(x), collapse = " | ")
      )
    }
    names(x)[hit]
  }

  species_col <- find_s3_column("species name")
  genome_col <- find_s3_column("genome type")
  novelty_col <- find_s3_column("novelty")
  pathogen_col <- find_s3_column("pathogen characteristics")
  kingdom_col <- find_s3_column("kingdom")

  # Remove any completely empty trailing rows before matching.
  keep <- !is.na(x[[species_col]]) & nzchar(trimws(as.character(x[[species_col]])))
  x <- x[keep, , drop = FALSE]
  x$match_name <- normalize_name(x[[species_col]])
  if (anyDuplicated(x$match_name)) {
    stop("Duplicated operational names in Table S3: ",
         paste(unique(x$match_name[duplicated(x$match_name)]), collapse = ", "))
  }
  kingdom <- tolower(strip_rank_prefix(x[[kingdom_col]]))
  genome <- tolower(as.character(x[[genome_col]]))
  x$Classification <- ifelse(
    kingdom == "bacteria", "Bacteria",
    ifelse(
      kingdom == "eukaryota", "Eukaryota",
      ifelse(grepl("dna", genome) & !grepl("ssrna-rt", genome), "DNA virus", "RNA virus")
    )
  )
  x$Known_pathogen <- x[[pathogen_col]] == "Known pathogen"
  x
}

read_organ_summary <- function(path, organ, s3) {
  sheet <- find_summary_sheet(path)
  dat <- read.xlsx(path, sheet = sheet, check.names = FALSE)
  if (ncol(dat) < 4) stop("Summary sheet has fewer than four columns: ", basename(path))

  abundance_col <- find_one_column(names(dat), "average log10|non zero avg log", path)
  count_col <- find_one_column(names(dat), "positive individuals|non zero count", path)

  # The first two columns are full operational name and display label.
  full_original <- trimws(as.character(dat[[1]]))
  display_original <- trimws(as.character(dat[[2]]))
  full_name <- normalize_name(full_original)
  display_name <- normalize_name(display_original)
  renamed_full <- full_name != full_original
  display_name[renamed_full & display_original == full_original] <- full_name[
    renamed_full & display_original == full_original
  ]

  idx <- match(full_name, s3$match_name)
  if (anyNA(idx)) {
    stop(
      "Unmatched Table S3 name(s) in ", basename(path), ": ",
      paste(full_original[is.na(idx)], collapse = ", ")
    )
  }

  out <- data.frame(
    Organ = organ,
    Input_row = seq_len(nrow(dat)) + 1L,
    Full_name = full_name,
    Display_name = display_name,
    Abundance = suppressWarnings(as.numeric(dat[[abundance_col]])),
    Positive_count = suppressWarnings(as.numeric(dat[[count_col]])),
    Classification = factor(
      s3$Classification[idx],
      levels = names(classification_colors)
    ),
    Known_pathogen = s3$Known_pathogen[idx],
    stringsAsFactors = FALSE
  )

  if (any(!is.finite(out$Abundance)) || any(!is.finite(out$Positive_count))) {
    stop("Non-numeric or missing summary values in ", basename(path))
  }
  if (any(out$Abundance < 0) || any(out$Positive_count < 0)) {
    stop("Negative abundance/count in ", basename(path))
  }
  if (anyDuplicated(out$Display_name)) {
    stop("Duplicated display names in ", basename(path), ": ",
         paste(unique(out$Display_name[duplicated(out$Display_name)]), collapse = ", "))
  }

  # Do not arrange or sort: this factor reverses only the plotting direction so
  # row 1 appears at the top while retaining the Excel summary order exactly.
  out$Display_factor <- factor(out$Display_name, levels = rev(out$Display_name))
  out
}

# ---------------------------- plot construction -----------------------------
make_organ_plot <- function(dat, organ) {
  count_top <- max(1, ceiling(max(dat$Positive_count, na.rm = TRUE)))
  abundance_pretty <- pretty(c(0, max(dat$Abundance, na.rm = TRUE)), n = 4)
  abundance_top <- max(abundance_pretty[is.finite(abundance_pretty)])
  if (!is.finite(abundance_top) || abundance_top <= 0) abundance_top <- 1
  scale_factor <- count_top / abundance_top
  dat$Abundance_scaled <- dat$Abundance * scale_factor

  ggplot(dat, aes(y = Display_factor)) +
    geom_col(
      aes(x = Positive_count, fill = Classification),
      width = bar_width,
      alpha = bar_alpha,
      color = NA
    ) +
    geom_path(
      aes(x = Abundance_scaled, group = 1),
      color = abundance_color,
      linewidth = line_width
    ) +
    geom_point(
      aes(x = Abundance_scaled),
      color = abundance_color,
      size = point_size,
      stroke = 0
    ) +
    # Draw every species label manually so known pathogens can be red while all
    # other labels remain black. Anchor at x = 0 (inside the scale) and extend
    # left with hjust; negative x values would be discarded by scale limits.
    geom_text(
      aes(
        x = 0,
        label = Display_name,
        color = Known_pathogen
      ),
      inherit.aes = TRUE,
      family = font_family,
      fontface = "italic",
      size = species_font_size / ggplot2::.pt,
      hjust = 1.05,
      vjust = 0.4,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = classification_colors, drop = FALSE) +
    scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = known_star_color),
      guide = "none"
    ) +
    scale_y_discrete(labels = NULL, expand = expansion(add = c(0.25, 0.25))) +
    scale_x_continuous(
      name = "Number of positive individuals",
      position = "top",
      limits = c(0, count_top),
      breaks = nice_integer_breaks(count_top),
      expand = expansion(mult = c(0, 0)),
      sec.axis = sec_axis(
        ~ . / scale_factor,
        name = expression("Average abundance " * log[10] * "(RPM+1)"),
        breaks = abundance_pretty[abundance_pretty >= 0 & abundance_pretty <= abundance_top]
      )
    ) +
    labs(title = organ, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_classic(base_family = font_family, base_size = axis_font_size) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        family = font_family, face = "bold", size = organ_title_size,
        hjust = 0, margin = margin(b = 2)
      ),
      axis.text.y = element_blank(),
      axis.text.x.top = element_text(
        family = font_family, size = axis_font_size - 0.5,
        color = "black", margin = margin(b = 1)
      ),
      axis.text.x.bottom = element_text(
        family = font_family, size = axis_font_size - 0.5,
        color = "black", margin = margin(t = 1)
      ),
      axis.title.x.top = element_text(
        family = font_family, size = axis_font_size,
        color = "black", margin = margin(b = 2)
      ),
      axis.title.x.bottom = element_text(
        family = font_family, size = axis_font_size,
        color = "black", margin = margin(t = 2)
      ),
      axis.ticks.y = element_blank(),
      axis.line.y = element_line(color = "black", linewidth = 0.35),
      axis.line.x = element_line(color = "black", linewidth = 0.35),
      axis.ticks.x = element_line(color = "black", linewidth = 0.3),
      plot.margin = margin(t = 5, r = 7, b = 5, l = 110)
    )
}

make_shared_legend <- function(horizontal = TRUE) {
  if (horizontal) {
    grobTree(
      textGrob("Data type", x = 0.01, y = 0.80, just = "left",
               gp = gpar(fontfamily = font_family, fontface = "bold", fontsize = 10)),
      rectGrob(x = 0.015, y = 0.44, width = 0.022, height = 0.22,
               just = "left", gp = gpar(fill = "black", col = NA)),
      textGrob("Number of positive individuals", x = 0.042, y = 0.44, just = "left",
               gp = gpar(fontfamily = font_family, fontsize = 9)),
      segmentsGrob(x0 = 0.235, x1 = 0.270, y0 = 0.44, y1 = 0.44,
                   gp = gpar(col = abundance_color, lwd = 1.2)),
      pointsGrob(x = 0.2525, y = 0.44, pch = 16, size = unit(2.0, "mm"),
                 gp = gpar(col = abundance_color, fill = abundance_color)),
      textGrob(expression("Average abundance " * log[10] * "(RPM+1)"),
               x = 0.278, y = 0.44, just = "left",
               gp = gpar(fontfamily = font_family, fontsize = 9)),

      textGrob("Classification", x = 0.42, y = 0.80, just = "left",
               gp = gpar(fontfamily = font_family, fontface = "bold", fontsize = 10)),
      rectGrob(x = c(0.425, 0.525, 0.620, 0.715), y = rep(0.44, 4),
               width = 0.018, height = 0.22, just = "left",
               gp = gpar(fill = unname(classification_colors), col = NA)),
      textGrob(names(classification_colors),
               x = c(0.448, 0.548, 0.643, 0.738), y = rep(0.44, 4), just = "left",
               gp = gpar(fontfamily = font_family, fontsize = 8.5)),

      textGrob("Known pathogens", x = 0.83, y = 0.80, just = "left",
               gp = gpar(fontfamily = font_family, fontface = "bold", fontsize = 10)),
      textGrob("Known pathogen species", x = 0.84, y = 0.44, just = "left",
               gp = gpar(fontfamily = font_family, fontface = "italic",
                         col = known_star_color, fontsize = 8.5))
    )
  } else {
    grobTree(
      textGrob("Data type", x = 0.05, y = 0.95, just = c("left", "top"),
               gp = gpar(fontfamily = font_family, fontface = "bold", fontsize = 11)),
      rectGrob(x = 0.06, y = 0.84, width = 0.08, height = 0.045,
               just = "left", gp = gpar(fill = "black", col = NA)),
      textGrob("Number of positive individuals", x = 0.17, y = 0.84, just = "left",
               gp = gpar(fontfamily = font_family, fontsize = 9)),
      segmentsGrob(x0 = 0.06, x1 = 0.14, y0 = 0.76, y1 = 0.76,
                   gp = gpar(col = abundance_color, lwd = 1.2)),
      pointsGrob(x = 0.10, y = 0.76, pch = 16, size = unit(2.0, "mm"),
                 gp = gpar(col = abundance_color, fill = abundance_color)),
      textGrob(expression("Average abundance " * log[10] * "(RPM+1)"),
               x = 0.17, y = 0.76, just = "left",
               gp = gpar(fontfamily = font_family, fontsize = 9)),
      textGrob("Classification", x = 0.05, y = 0.66, just = "left",
               gp = gpar(fontfamily = font_family, fontface = "bold", fontsize = 11)),
      rectGrob(x = rep(0.06, 4), y = c(0.57, 0.49, 0.41, 0.33),
               width = 0.08, height = 0.045, just = "left",
               gp = gpar(fill = unname(classification_colors), col = NA)),
      textGrob(names(classification_colors), x = rep(0.17, 4),
               y = c(0.57, 0.49, 0.41, 0.33), just = "left",
               gp = gpar(fontfamily = font_family, fontsize = 9)),
      textGrob("Known pathogens", x = 0.05, y = 0.22, just = "left",
               gp = gpar(fontfamily = font_family, fontface = "bold", fontsize = 11)),
      textGrob("Known pathogen species", x = 0.06, y = 0.12, just = "left",
               gp = gpar(fontfamily = font_family, fontface = "italic",
                         col = known_star_color, fontsize = 9))
    )
  }
}

save_grob_pdf <- function(grob, filename, width, height) {
  grDevices::cairo_pdf(filename, width = width, height = height, family = font_family)
  grid.newpage()
  grid.draw(grob)
  dev.off()
}

# ---------------------------- read, plot, save -------------------------------
s3 <- read_s3(s3_file)
organ_data <- setNames(vector("list", length(organ_order)), organ_order)
plots <- setNames(vector("list", length(organ_order)), organ_order)

for (organ in organ_order) {
  organ_data[[organ]] <- read_organ_summary(organ_files[[organ]], organ, s3)
  plots[[organ]] <- make_organ_plot(organ_data[[organ]], organ)
}

# Convert to gtables and enforce identical widths before arranging 2 rows × 4 columns.
plot_gtables <- lapply(plots[organ_order], ggplotGrob)
common_widths <- Reduce(grid::unit.pmax, lapply(plot_gtables, `[[`, "widths"))
plot_gtables <- lapply(plot_gtables, function(g) {
  g$widths <- common_widths
  g
})

row_max_taxa <- c(
  max(vapply(organ_data[c("Brain", "Gill", "Intestine", "Kidney")], nrow, integer(1))),
  max(vapply(organ_data[c("Liver", "Muscle", "Skin", "Spleen")], nrow, integer(1)))
)
combined_grid <- arrangeGrob(
  grobs = plot_gtables,
  ncol = 4,
  nrow = 2,
  widths = rep(1, 4),
  heights = row_max_taxa
)
combined <- arrangeGrob(
  make_shared_legend(horizontal = TRUE),
  combined_grid,
  ncol = 1,
  heights = c(0.065, 0.935)
)
save_grob_pdf(
  combined,
  file.path(script_dir, "figure3_combined_2rows_4cols.pdf"),
  width = combined_width_in,
  height = combined_height_in
)

# Machine-readable audit of the exact plotted rows, in unchanged input order.
plot_audit <- do.call(rbind, lapply(organ_data, function(x) {
  x[, c(
    "Organ", "Input_row", "Full_name", "Display_name", "Abundance",
    "Positive_count", "Classification", "Known_pathogen"
  )]
}))
write.table(
  plot_audit,
  file = file.path(script_dir, "figure3_plot_mapping.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8"
)
