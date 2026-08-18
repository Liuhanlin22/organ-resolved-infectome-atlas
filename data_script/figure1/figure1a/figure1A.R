#!/usr/bin/env Rscript

# Figure 1A — COX1 host phylogeny
#
# Input annotation table (tab-separated; header required):
# original_id, display_species, family
#
# Run (edit the two paths below, or provide them as command-line arguments):
# Rscript figure1A.R cox1_rooted.treefile.nwk cox1_tree_annotation.tsv figure1A
#
# The script exports a publication-quality PDF and PNG.  Tips with an ID that
# contains "YS-" are the study sequences and receive species-specific dots and
# matching pale backgrounds.

options(stringsAsFactors = FALSE)

required_cran <- c("ape", "ggplot2", "dplyr", "readr", "tibble", "grid")
required_bioc <- c("ggtree", "treeio")

install_if_missing <- function() {
  cran_missing <- required_cran[!vapply(required_cran, requireNamespace,
                                        logical(1), quietly = TRUE)]
  if (length(cran_missing) > 0) {
    install.packages(cran_missing, repos = "https://cloud.r-project.org")
  }
  bioc_missing <- required_bioc[!vapply(required_bioc, requireNamespace,
                                        logical(1), quietly = TRUE)]
  if (length(bioc_missing) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
  }
}
install_if_missing()

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

# Default paths: modify only if you do not use command-line arguments.
tree_file_default <- "cox1_rooted.treefile.nwk"
annotation_file_default <- "cox1_tree_annotation.tsv"
output_prefix_default <- "figure1A"

tree_file <- if (length(args) >= 1) args[1] else tree_file_default
annotation_file <- if (length(args) >= 2) args[2] else annotation_file_default
output_prefix <- if (length(args) >= 3) args[3] else output_prefix_default

if (!file.exists(tree_file)) stop("Tree file not found: ", tree_file)
if (!file.exists(annotation_file)) stop("Annotation file not found: ", annotation_file)

# Appearance parameters (adjust here when needed)
tip_size <- 2.8
family_size <- 3.20
tree_panel_width <- 0.2

tree <- ape::read.tree(tree_file)
ann <- readr::read_tsv(annotation_file, show_col_types = FALSE,
                       progress = FALSE, name_repair = "minimal") %>%
  select(original_id, display_species, family) %>%
  mutate(
    original_id = trimws(original_id),
    display_species = trimws(display_species),
    family = trimws(family),
    sequence_source = if_else(grepl("YS-", original_id, fixed = TRUE),
                              "Study sequence", "Reference sequence")
  )

if (anyDuplicated(ann$original_id)) {
  duplicated_ids <- unique(ann$original_id[duplicated(ann$original_id)])
  stop("Duplicate original_id values in annotation table: ",
       paste(duplicated_ids, collapse = "; "))
}

missing_annotation <- setdiff(tree$tip.label, ann$original_id)
extra_annotation <- setdiff(ann$original_id, tree$tip.label)
if (length(missing_annotation) > 0) {
  stop("The following tree tip(s) have no annotation: \n",
       paste(missing_annotation, collapse = "\n"))
}
if (length(extra_annotation) > 0) {
  warning("Annotation rows not present in tree (ignored): \n",
          paste(extra_annotation, collapse = "\n"))
}
if (any(ann$display_species == "" | is.na(ann$display_species)) ||
    any(ann$family == "" | is.na(ann$family))) {
  stop("display_species and family must be filled for every sequence.")
}

# Attach metadata to ggtree data.  The internal-node numeric labels are retained
# as ultrafast-bootstrap values from IQ-TREE.
p <- ggtree(tree, size = 0.15, colour = "#222222")
p$data <- p$data %>%
  left_join(ann, by = c("label" = "original_id"))

tip_data <- p$data %>%
  filter(isTip) %>%
  select(node, label, x, y, display_species, family, sequence_source)

if (any(is.na(tip_data$display_species)) || any(is.na(tip_data$family))) {
  stop("Metadata join failed for at least one tip; check original_id exactly matches tree labels.")
}

# A single right-side bracket is meaningful only when a family's tips are one
# continuous block in the displayed tree.  Warn about non-contiguous families;
# their bracket still shows the complete y-range, so these should be checked.
family_ranges <- tip_data %>%
  arrange(y) %>%
  group_by(family) %>%
  summarise(
    ymin = min(y), ymax = max(y), n_tips = n(),
    contiguous = (max(y) - min(y) + 1 == n()),
    .groups = "drop"
  ) %>%
  arrange(desc(ymax))

# Each study species receives its own colour.  Its dot uses the saturated
# colour and its label background uses a light version of that same colour.
# Reference-sequence labels remain black, keeping the study samples prominent.
study_species <- tip_data %>%
  filter(sequence_source == "Study sequence") %>%
  distinct(display_species) %>%
  arrange(display_species) %>%
  pull(display_species)
species_palette <- setNames(
  grDevices::hcl(
    h = seq(10, 370, length.out = length(study_species) + 1)[-1],
    c = 72, l = 48
  ),
  study_species
)
tip_data <- tip_data %>%
  mutate(
    species_colour = unname(species_palette[display_species]),
    species_fill = grDevices::adjustcolor(species_colour, alpha.f = 0.1)
  )

non_contiguous <- family_ranges %>% filter(!contiguous)
if (nrow(non_contiguous) > 0) {
  warning(
    "Some family labels are not contiguous in this gene tree: ",
    paste(non_contiguous$family, collapse = "; "),
    ". Their right-side bracket spans intervening taxa; inspect these manually."
  )
}

# This fixed-width canvas makes `tree_panel_width` a true visual-width control.
# Merely multiplying all branch lengths does not work: ggplot rescales the
# complete x-axis and the plotted tree keeps the same apparent width.
x_tree_max <- max(tip_data$x)
x_range <- diff(range(p$data$x))
if (x_range == 0) x_range <- 1
canvas_xmax <- x_tree_max / tree_panel_width
tip_x <- x_tree_max + 0.025 * x_range
dot_x <- tip_x - 0.020 * x_range
family_x <- 0.31 * canvas_xmax
bracket_tick <- 0.025 * x_range

family_ranges <- family_ranges %>%
  mutate(y_mid = (ymin + ymax) / 2)

# Reference sequences are plain italic text. Study sequences are labels with a
# background box, which is more legible than colouring the entire text line.
p <- p +
  geom_segment(
  data = tip_data,
  aes(x = x, y = y, yend = y),
  xend = dot_x,
  inherit.aes = FALSE,
  linetype = "dashed",
  linewidth = 0.25,
  colour = "#777777"
) +
  geom_point(
    data = tip_data %>% filter(sequence_source == "Study sequence"),
    aes(x = dot_x, y = y, colour = display_species),
    inherit.aes = FALSE, shape = 16, size = 2.30
  ) +
  geom_text(
    data = tip_data %>% filter(sequence_source == "Reference sequence"),
    aes(x = tip_x, y = y, label = display_species),
    inherit.aes = FALSE, hjust = 0, fontface = "italic", size = tip_size
  ) +
  geom_label(
    data = tip_data %>% filter(sequence_source == "Study sequence"),
    aes(x = tip_x, y = y, label = display_species, fill = species_fill),
    inherit.aes = FALSE, hjust = 0, fontface = "italic", size = tip_size,
    colour = "#111111", label.size = 0,
    label.padding = unit(c(0.04, 0.10, 0.04, 0.10), "lines")
  ) +
  geom_segment(
    data = family_ranges,
    aes(x = family_x, xend = family_x, y = ymin - 0.35, yend = ymax + 0.35),
    inherit.aes = FALSE, colour = "black", linewidth = 0.2
  ) +
  geom_segment(
    data = family_ranges,
    aes(x = family_x - bracket_tick, xend = family_x, y = ymin - 0.35, yend = ymin - 0.35),
    inherit.aes = FALSE, colour = "black", linewidth = 0.2
  ) +
  geom_segment(
    data = family_ranges,
    aes(x = family_x - bracket_tick, xend = family_x, y = ymax + 0.35, yend = ymax + 0.35),
    inherit.aes = FALSE, colour = "black", linewidth = 0.2
  ) +
  geom_text(
    data = family_ranges,
    aes(x = family_x + 0.020 * x_range, y = y_mid, label = family),
    inherit.aes = FALSE, hjust = 0, fontface = "plain", size = family_size,
    colour = "black"
  )

# Put a manual scale bar in the upper-left using the original branch-length unit.
scale_x0 <- 0.020 * x_range
scale_width <- 0.05
scale_y <- max(tip_data$y) - 4.0
p <- p +
  annotate("text", x = scale_x0, y = scale_y + 1.55, label = "COX1",
           hjust = 0, size = 3.8) +
  annotate("segment", x = scale_x0, xend = scale_x0 + scale_width,
           y = scale_y, yend = scale_y, linewidth = 0.2) +
  annotate("text", x = scale_x0 + scale_width / 2, y = scale_y - 1.15,
           label = "0.05", hjust = 0.5, size = 3.0)

# The fixed x-limit makes tree_panel_width affect the actual rendered width.
p <- p +
  xlim(0, canvas_xmax) +
  scale_colour_manual(values = species_palette, guide = "none") +
  scale_fill_identity() +
  theme_tree() +
  theme(
    text = element_text(family = "Arial"),
    plot.margin = margin(8, 20, 8, 8)
  )

pdf_file <- paste0(output_prefix, ".pdf")
png_file <- paste0(output_prefix, ".png")

# The dynamic canvas is compact while retaining an enlarged, readable species
# font.  Adjust 0.09 downward slightly if the final export still has excess gap.
fig_height <- max(14, nrow(tip_data) * 0.09 + 1.5)
ggsave(pdf_file, p, width = 18, height = fig_height, units = "in",
       limitsize = FALSE, device = grDevices::pdf)
ggsave(png_file, p, width = 18, height = fig_height, units = "in", dpi = 600,
       limitsize = FALSE, bg = "white")

