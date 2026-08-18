#!/usr/bin/env Rscript

# Figure 4B — organ distribution of selected pathogens
# Input:  figure4A.xlsx
# Output: Figure4B.pdf, Figure4B.png and verification tables

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

suppressWarnings(
  suppressPackageStartupMessages({
    library(openxlsx)
    library(ggplot2)
    library(gridExtra)
  })
)

# ---------------------------- adjustable parameters -------------------------
font_family <- "Arial"
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
}

# Species-title typography. Change these four values to modify every line of
# every species title at once. Each title below is ONE character string; "\n"
# only inserts a line break and does not create separate R text grobs.
species_title_font_family <- font_family
species_title_font_size <- 9.0
species_title_font_face <- "italic"
species_title_lineheight <- 0.90

figure_width_in <- 11.0
figure_height_in <- 8.2
bubble_size <- 8.4
bubble_stroke <- 0.30
bubble_column_width <- 0.75
first_column_width <- 0.94
count_column_width <- 1.00
panel_gap_width <- 0.15

organ_order <- c(
  "Brain", "Gill", "Intestine", "Liver",
  "Kidney", "Spleen", "Muscles", "Skin"
)

organ_colors <- c(
  "Brain"     = "#FEC6DC",
  "Gill"      = "#FD8701",
  "Intestine" = "#AB696B",
  "Liver"     = "#A076E3",
  "Kidney"    = "#D6BFEA",
  "Spleen"    = "#0095C5",
  "Muscles"   = "#ABE7F2",
  "Skin"      = "#B5E3C6"
)

# Panel order: first 8 are the upper row, last 7 are the lower row.
# Species titles are drawn above their corresponding bubble panels.
selected_species <- c(
  "Bombay duck fish bornavirus",
  "Odontamblyopus rubicundus hepacivirus",
  "Johnius belangerii hepacivirus",
  "Nanhai reo-like virus 7",
  "Nanhai retrovirus 1",
  "Nanhai circo-like virus 2",
  "Glugea plecoglossi",
  "Toxocara sp. CTSZ",
  "Telatrygon zugei hepacivirus",
  "Vibrio ponticus",
  "Aliivibrio fischeri",
  "Photobacterium sp.",
  "Vibrio sp. Y",
  "Asian seabass Nervous Necrosis Virus",
  "Marble sleepy goby iridovirus"
)

# Compact display labels matching the reference layout.
# IMPORTANT: keep each value as one quoted string and use \n for line breaks.
# Example: "Bombay\nduck fish\nbornavirus" is one editable title in the script.
species_display_labels <- c(
  "Bombay duck fish bornavirus" = "Bombay\nduck fish\nbornavirus",
  "Odontamblyopus rubicundus hepacivirus" = "Odont.\nrubicundus\nhepacivirus",
  "Johnius belangerii hepacivirus" = "Johnius\nbelangerii\nhepacivirus",
  "Nanhai reo-like virus 7" = "Nanhai\nreo-like\nvirus 7",
  "Nanhai retrovirus 1" = "Nanhai\nretrovirus\n1",
  "Nanhai circo-like virus 2" = "Nanhai\ncirco-like\nvirus 2",
  "Glugea plecoglossi" = "Glugea\nplecoglossi",
  "Toxocara sp. CTSZ" = "Toxocara\nsp. CTSZ",
  "Telatrygon zugei hepacivirus" = "Telatrygon\nzugei\nhepacivirus",
  "Vibrio ponticus" = "Vibrio\nponticus",
  "Aliivibrio fischeri" = "Aliivibrio\nfischeri",
  "Photobacterium sp." = "Photo\nbacterium\nsp.",
  "Vibrio sp. Y" = "Vibrio\nsp. Y",
  "Asian seabass Nervous Necrosis Virus" = "Asian seabass\nNervous\nNecrosis virus",
  "Marble sleepy goby iridovirus" = "Marble\nsleepy goby\niridovirus"
)

# Red-starred panels in the supplied reference figure.
starred_species <- c(
  "Bombay duck fish bornavirus",
  "Glugea plecoglossi",
  "Vibrio ponticus",
  "Aliivibrio fischeri",
  "Asian seabass Nervous Necrosis Virus",
  "Marble sleepy goby iridovirus"
)

bar_fill <- "#D9D2E6"
bar_left_limit <- 10
bar_right_min <- 38
bar_right_max <- 42
bar_right_start <- 11.6

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
input_file <- file.path(script_dir, "figure4A.xlsx")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

normalize_header <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", " ", x)
  trimws(gsub("[[:space:]]+", " ", x))
}

# ---------------------------- read and validate ------------------------------
dat <- read.xlsx(input_file, sheet = 1, check.names = FALSE)
original_names <- names(dat)
header_keys <- normalize_header(original_names)

find_column <- function(key) {
  hit <- which(header_keys == normalize_header(key))
  if (length(hit) != 1) {
    stop(
      "Column not found or duplicated: ", key,
      "\nDetected columns: ", paste(original_names, collapse = " | ")
    )
  }
  original_names[hit]
}

pathogen_col <- find_column("Pathogens")
organ_cols <- setNames(vapply(organ_order, find_column, character(1)), organ_order)

dat$Pathogens_for_plot <- trimws(as.character(dat[[pathogen_col]]))
if (anyNA(dat$Pathogens_for_plot) || any(!nzchar(dat$Pathogens_for_plot))) {
  stop("The Pathogens column contains empty values.")
}
if (anyDuplicated(dat$Pathogens_for_plot)) {
  duplicates <- unique(dat$Pathogens_for_plot[duplicated(dat$Pathogens_for_plot)])
  stop("Duplicated pathogen names: ", paste(duplicates, collapse = ", "))
}

for (organ in organ_order) {
  values <- suppressWarnings(as.numeric(dat[[organ_cols[[organ]]]]))
  if (any(!is.finite(values))) {
    stop("Non-numeric or missing values detected in organ column: ", organ)
  }
  if (any(values < 0)) stop("Negative values detected in organ column: ", organ)
  dat[[organ]] <- values
}

missing_species <- setdiff(selected_species, dat$Pathogens_for_plot)
if (length(missing_species) > 0) {
  stop("Selected pathogen(s) not found: ", paste(missing_species, collapse = ", "))
}

selected <- dat[match(selected_species, dat$Pathogens_for_plot), , drop = FALSE]

# ---------------------------- bubble-plot data -------------------------------
bubble_data <- do.call(rbind, lapply(seq_along(selected_species), function(i) {
  values <- as.numeric(unlist(selected[i, organ_order, drop = FALSE], use.names = FALSE))
  data.frame(
    Panel = if (i <= 8) paste0("Top_", i) else paste0("Bottom_", i - 8),
    Panel_index = i,
    Pathogens = selected_species[i],
    Organ = organ_order,
    RPM = values,
    Log10_RPM_plus_1 = log10(values + 1),
    stringsAsFactors = FALSE
  )
}))
bubble_data$Organ <- factor(bubble_data$Organ, levels = rev(organ_order))

panel_map <- unique(bubble_data[, c("Panel", "Panel_index", "Pathogens")])
panel_map$Display_label <- unname(species_display_labels[panel_map$Pathogens])
panel_map$Red_star <- panel_map$Pathogens %in% starred_species
write.table(
  panel_map,
  file = file.path(script_dir, "Figure4B_panel_species_mapping.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8"
)
write.table(
  bubble_data,
  file = file.path(script_dir, "Figure4B_selected_species_values.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8"
)

make_bubble_plot <- function(pathogen, show_y_labels = FALSE) {
  d <- bubble_data[bubble_data$Pathogens == pathogen, , drop = FALSE]
  x_upper <- max(1, ceiling(max(d$Log10_RPM_plus_1)))
  x_breaks <- seq(0, x_upper, by = 1)
  x_minor_breaks <- if (x_upper <= 2) {
    seq(0.5, x_upper - 0.5, by = 1)
  } else {
    NULL
  }
  band_data <- data.frame(
    xmin = seq(0, x_upper - 1, by = 1),
    xmax = seq(0, x_upper - 1, by = 1) + 0.5
  )

  y_text <- if (show_y_labels) {
    element_text(
      family = font_family, size = 10.5, color = "black",
      margin = margin(r = 4)
    )
  } else {
    element_blank()
  }

  ggplot(d, aes(x = Log10_RPM_plus_1, y = Organ)) +
    geom_rect(
      data = band_data,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "#D9D9D9",
      color = NA,
      alpha = 0.78
    ) +
    geom_point(
      aes(fill = Organ),
      shape = 21,
      size = bubble_size,
      stroke = bubble_stroke,
      color = "black"
    ) +
    scale_fill_manual(values = organ_colors, guide = "none") +
    scale_x_continuous(
      limits = c(0, x_upper),
      breaks = x_breaks,
      minor_breaks = x_minor_breaks,
      labels = rep("", length(x_breaks)),
      sec.axis = dup_axis(
        breaks = x_breaks,
        labels = x_breaks,
        name = NULL
      ),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = expression(log[10](RPM + 1)),
      y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_family = font_family, base_size = 8.5) +
    theme(
      plot.title = element_blank(),
      panel.border = element_rect(color = "#666666", fill = NA, linewidth = 0.4),
      panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.35),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_line(
        color = "#E7E7E7", linewidth = 0.30
      ),
      axis.text.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank(),
      axis.text.x.top = element_text(
        family = font_family, size = 8.5, color = "black",
        margin = margin(b = 2)
      ),
      axis.ticks.x.top = element_line(color = "black", linewidth = 0.35),
      axis.title.x = element_text(
        family = font_family, size = 8.2, color = "black",
        margin = margin(t = 3)
      ),
      axis.text.y = y_text,
      axis.ticks.y = element_line(color = "#BFBFBF", linewidth = 0.3),
      plot.margin = margin(t = 0, r = 2.0, b = 1.5, l = 2.0)
    )
}

make_bubble_component <- function(pathogen, show_y_labels = FALSE) {
  title_label <- unname(species_display_labels[pathogen])
  if (is.na(title_label) || !nzchar(title_label)) title_label <- pathogen
  title_line_count <- length(strsplit(title_label, "\n", fixed = TRUE)[[1]])

  title_text <- grid::textGrob(
    label = title_label,
    x = grid::unit(if (pathogen %in% starred_species) 0.55 else 0.50, "npc"),
    y = grid::unit(0.50, "npc"),
    just = "center",
    gp = grid::gpar(
      fontfamily = species_title_font_family,
      fontsize = species_title_font_size,
      fontface = species_title_font_face,
      lineheight = species_title_lineheight,
      col = "black"
    )
  )

  if (pathogen %in% starred_species) {
    star_x <- if (identical(pathogen, selected_species[1])) 0.32 else 0.20
    star_y <- if (title_line_count >= 3) 0.76 else 0.66
    title_grob <- grid::grobTree(
      title_text,
      grid::textGrob(
        "★",
        x = grid::unit(star_x, "npc"),
        y = grid::unit(star_y, "npc"),
        just = "left",
        gp = grid::gpar(
          fontfamily = font_family,
          fontsize = 10.5,
          fontface = "plain",
          col = "red"
        )
      )
    )
  } else {
    title_grob <- title_text
  }

  arrangeGrob(
    title_grob,
    make_bubble_plot(pathogen, show_y_labels = show_y_labels),
    ncol = 1,
    heights = c(0.24, 1)
  )
}

# ---------------------------- infected-organ count --------------------------
# For every pathogen, count the number of organs with RPM > 0. Then count how
# many pathogen species occur at each breadth (1 to 8 infected organs).
positive_matrix <- as.matrix(dat[, organ_order, drop = FALSE]) > 0
infected_n <- rowSums(positive_matrix)
if (any(infected_n < 1 | infected_n > 8)) {
  stop("Every pathogen must occur in 1 to 8 organs.")
}

infected_count <- as.data.frame(table(factor(infected_n, levels = 1:8)))
names(infected_count) <- c("Number_of_infected_organs", "Pathogen_species_count")
infected_count$Number_of_infected_organs <- as.integer(
  as.character(infected_count$Number_of_infected_organs)
)
infected_count$Organ_factor <- factor(
  infected_count$Number_of_infected_organs,
  levels = rev(1:8)
)

write.table(
  infected_count[, c("Number_of_infected_organs", "Pathogen_species_count")],
  file = file.path(script_dir, "Figure4B_infected_organ_count.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8"
)

max_count <- max(infected_count$Pathogen_species_count)

# Draw the broken axis inside one ggplot panel so that the two bar segments
# share exactly the same y coordinates. This prevents vertical misalignment.
count_left <- infected_count
count_left$xmin <- 0
count_left$xmax <- pmin(count_left$Pathogen_species_count, bar_left_limit)
count_left$ymin <- count_left$Number_of_infected_organs - 0.36
count_left$ymax <- count_left$Number_of_infected_organs + 0.36

count_right <- infected_count[
  infected_count$Pathogen_species_count >= bar_right_min,
  , drop = FALSE
]
count_right$xmin <- bar_right_start
count_right$xmax <- bar_right_start +
  (count_right$Pathogen_species_count - bar_right_min)
count_right$ymin <- count_right$Number_of_infected_organs - 0.36
count_right$ymax <- count_right$Number_of_infected_organs + 0.36

count_x_breaks <- c(
  seq(0, bar_left_limit, by = 2),
  bar_right_start + seq(0, bar_right_max - bar_right_min, by = 2)
)
count_x_labels <- c(
  seq(0, bar_left_limit, by = 2),
  seq(bar_right_min, bar_right_max, by = 2)
)
count_x_max <- bar_right_start + (bar_right_max - bar_right_min)

count_plot <- ggplot() +
  geom_rect(
    data = count_left,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = bar_fill,
    color = NA
  ) +
  geom_rect(
    data = count_right,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = bar_fill,
    color = NA
  ) +
  annotate(
    "segment",
    x = c(10.35, 10.78), xend = c(10.62, 11.05),
    y = c(0.52, 0.52), yend = c(0.76, 0.76),
    linewidth = 0.45, color = "black"
  ) +
  scale_x_continuous(
    position = "top",
    limits = c(0, count_x_max),
    breaks = count_x_breaks,
    labels = count_x_labels,
    sec.axis = dup_axis(
      breaks = count_x_breaks,
      labels = count_x_labels,
      name = NULL
    ),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_reverse(
    limits = c(8.5, 0.5),
    breaks = 1:8,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Number of infected organs"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = font_family, base_size = 8) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(family = font_family, color = "black", size = 7.2),
    axis.title.y = element_text(
      family = font_family, color = "black", size = 7.7,
      margin = margin(r = 1)
    ),
    axis.text.x.bottom = element_text(
      family = font_family, color = "black", size = 7.2
    ),
    axis.ticks.x.bottom = element_line(color = "black", linewidth = 0.3),
    panel.grid.major.x = element_line(
      color = "#E2E2E2", linewidth = 0.3, linetype = "dashed"
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 0, r = 4, b = 1.5, l = 4)
  )

make_count_component <- function() {
  count_title <- grid::textGrob(
    "Number of\npathogen\nspecies",
    x = grid::unit(0.50, "npc"),
    y = grid::unit(0.50, "npc"),
    just = "center",
    gp = grid::gpar(
      fontfamily = font_family,
      fontsize = 9.2,
      fontface = "plain",
      lineheight = 0.90,
      col = "black"
    )
  )

  arrangeGrob(
    count_title,
    count_plot,
    ncol = 1,
    heights = c(0.24, 1)
  )
}

# ---------------------------- assemble and export ---------------------------
interleave_grobs <- function(grobs) {
  result <- list()
  for (i in seq_along(grobs)) {
    result[[length(result) + 1]] <- grobs[[i]]
    if (i < length(grobs)) {
      result[[length(result) + 1]] <- grid::nullGrob()
    }
  }
  result
}

interleave_widths <- function(component_widths, gap_width) {
  result <- numeric(0)
  for (i in seq_along(component_widths)) {
    result <- c(result, component_widths[i])
    if (i < length(component_widths)) result <- c(result, gap_width)
  }
  result
}

build_final_grob <- function() {
  # Build grobs only after the output device is open. This allows Cairo/PNG
  # to resolve Arial directly and avoids PostScript font-database warnings.
  top_plots <- lapply(seq_len(8), function(i) {
    make_bubble_component(selected_species[i], show_y_labels = i == 1)
  })
  bottom_plots <- lapply(9:15, function(i) {
    make_bubble_component(selected_species[i], show_y_labels = i == 9)
  })

  top_components <- top_plots
  bottom_components <- c(bottom_plots, list(make_count_component()))
  top_component_widths <- c(first_column_width, rep(bubble_column_width, 7))
  bottom_component_widths <- c(
    first_column_width,
    rep(bubble_column_width, 6),
    count_column_width
  )
  top_right_spacer <- count_column_width - bubble_column_width

  top_grob <- arrangeGrob(
    grobs = c(interleave_grobs(top_components), list(grid::nullGrob())),
    ncol = 16,
    widths = c(
      interleave_widths(top_component_widths, panel_gap_width),
      top_right_spacer
    )
  )
  bottom_grob <- arrangeGrob(
    grobs = interleave_grobs(bottom_components),
    ncol = 15,
    widths = interleave_widths(bottom_component_widths, panel_gap_width)
  )

  arrangeGrob(
    top_grob,
    bottom_grob,
    ncol = 1,
    heights = c(1, 1)
  )
}

draw_figure <- function() {
  final_grob <- build_final_grob()
  grid::grid.draw(final_grob)
  grid::grid.text(
    "b",
    x = grid::unit(0.006, "npc"),
    y = grid::unit(0.992, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(
      fontfamily = font_family,
      fontsize = 20,
      fontface = "bold",
      col = "black"
    )
  )
}

pdf_file <- file.path(script_dir, "Figure4B.pdf")
grDevices::cairo_pdf(
  pdf_file,
  width = figure_width_in,
  height = figure_height_in,
  family = font_family,
  bg = "white"
)
draw_figure()
grDevices::dev.off()

png_file <- file.path(script_dir, "Figure4B.png")
grDevices::png(
  png_file,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  res = 600,
  bg = "white"
)
draw_figure()
grDevices::dev.off()
