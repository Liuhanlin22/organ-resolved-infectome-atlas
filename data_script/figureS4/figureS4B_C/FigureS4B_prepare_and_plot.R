# Figure S4B: pathogen counts by broad category and tissue
#
# The count is the number of classified pathogen taxa detected in at least one
# of the 112 complete libraries for each tissue. The four panels are RNA, DNA,
# BACTERIA and EUKARYOTA, in the tissue order used by the supplied reference.

required_packages <- c("ggplot2", "readxl", "dplyr", "tidyr", "openxlsx")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
})

script_arguments <- commandArgs(trailingOnly = FALSE)
script_file_argument <- grep("^--file=", script_arguments, value = TRUE)
script_dir <- if (length(script_file_argument) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_argument[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}

output_dir <- script_dir
abundance_file <- file.path(output_dir, "Abundance_matrix_with_libs_complete8_896.xlsx")
metadata_file <- file.path(output_dir, "Sample_metadata_complete8_896.xlsx")
annotation_file <- file.path(script_dir, "..", "..", "figure2", "figure2b", "row_heatmap_annot.txt")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tissue_order <- c("Intestine", "Spleen", "Kidney", "Skin", "Gill", "Liver", "Muscles", "Brain")
category_order <- c("RNA", "DNA", "BACTERIA", "EUKARYOTA")
source_to_category <- c(
  "RNA virus" = "RNA",
  "DNA virus" = "DNA",
  "Bacteria" = "BACTERIA",
  "Eukaryota" = "EUKARYOTA"
)
pathogen_aliases <- c(
  "Cystoisospora_sp._DSLT-2" = "Cystoisospora sp. DSLT-2",
  "Photobacterium sp. QLXSD" = "Photobacterium damselae QLXSD"
)
category_fill <- c(
  RNA = "#F2DCDB",
  DNA = "#B8CAE0",
  BACTERIA = "#C8B7D9",
  EUKARYOTA = "#CFE0B6"
)
category_border <- c(
  RNA = "#B45A19",
  DNA = "#5C91BD",
  BACTERIA = "#72519A",
  EUKARYOTA = "#78953F"
)

abundance <- read_excel(abundance_file, sheet = "abundance")
names(abundance)[1] <- "Pathogen"
abundance$Pathogen <- trimws(as.character(abundance$Pathogen))
sample_columns <- names(abundance)[-1]
if (nrow(abundance) != 75L || length(sample_columns) != 896L) {
  stop("Expected 75 pathogens and 896 libraries in the filtered abundance matrix.")
}
abundance_numeric <- abundance %>%
  select(all_of(sample_columns)) %>%
  mutate(across(everything(), as.numeric))
if (anyNA(abundance_numeric) || any(as.matrix(abundance_numeric) < 0)) {
  stop("The abundance matrix contains missing, non-numeric, or negative values.")
}
rpm <- as.matrix(abundance_numeric)
rownames(rpm) <- abundance$Pathogen

metadata <- read_excel(metadata_file, sheet = "main") %>%
  transmute(
    lib_index = trimws(as.character(lib_index)),
    tissue_types = trimws(as.character(tissue_types))
  )
if (nrow(metadata) != 896L || !setequal(metadata$lib_index, sample_columns)) {
  stop("The metadata does not match the 896 abundance libraries.")
}
if (any(table(metadata$tissue_types)[tissue_order] != 112L)) {
  stop("Each tissue must contain exactly 112 libraries.")
}

classification <- read.delim(annotation_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!"Pathogen taxa" %in% names(classification) && "Taxa" %in% names(classification)) {
  names(classification)[names(classification) == "Taxa"] <- "Pathogen taxa"
}
classification$Pathogen <- trimws(as.character(classification$Pathogen))
classification$Pathogen <- ifelse(
  classification$Pathogen %in% names(pathogen_aliases),
  unname(pathogen_aliases[classification$Pathogen]),
  classification$Pathogen
)
classification$Category <- unname(source_to_category[as.character(classification$Classification)])
if (nrow(classification) != 75L || anyNA(classification$Category)) {
  stop("The classification table must contain 75 rows with known categories.")
}
if (!setequal(classification$Pathogen, rownames(rpm))) {
  stop("Pathogen names do not match after applying aliases.")
}

detected_long <- bind_rows(lapply(seq_len(nrow(classification)), function(i) {
  pathogen <- classification$Pathogen[i]
  tissue_rows <- lapply(tissue_order, function(tissue) {
    libraries <- metadata$lib_index[metadata$tissue_types == tissue]
    positive_libraries <- sum(rpm[pathogen, libraries] > 0)
    tibble(
      Pathogen = pathogen,
      `Pathogen taxa` = classification$`Pathogen taxa`[i],
      Classification = classification$Classification[i],
      Category = classification$Category[i],
      Tissue = tissue,
      Positive_libraries = positive_libraries,
      Detected = positive_libraries > 0
    )
  })
  bind_rows(tissue_rows)
}))

counts_long <- detected_long %>%
  group_by(Category, Tissue) %>%
  summarise(
    Pathogen_count = sum(Detected),
    Definition = "Number of classified pathogens with RPM > 0 in at least one of 112 libraries",
    .groups = "drop"
  ) %>%
  mutate(
    Category = factor(Category, levels = category_order),
    Tissue = factor(Tissue, levels = tissue_order)
  ) %>%
  arrange(Category, Tissue)

counts_wide <- counts_long %>%
  select(Category, Tissue, Pathogen_count) %>%
  pivot_wider(names_from = Tissue, values_from = Pathogen_count) %>%
  arrange(Category)

write.csv(detected_long, file.path(output_dir, "FigureS4B_pathogen_detection_long.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(counts_long, file.path(output_dir, "FigureS4B_counts_long.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(counts_wide, file.path(output_dir, "FigureS4B_counts_wide.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.table(
  classification[, c("Pathogen", "Pathogen taxa", "Classification", "Category")],
  file.path(output_dir, "FigureS4B_pathogen_classification.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8"
)

plot_data <- counts_long %>%
  mutate(
    Category = factor(Category, levels = category_order),
    Tissue = factor(Tissue, levels = tissue_order),
    x = match(Tissue, tissue_order),
    panel_max = c(RNA = 20, DNA = 3, BACTERIA = 10, EUKARYOTA = 12)[as.character(Category)]
  )

# The colored rectangles provide the panel borders while allowing each panel
# to retain its own y-axis scale.
panel_data <- plot_data %>%
  distinct(Category, panel_max) %>%
  mutate(xmin = 0.45, xmax = 8.55, ymin = 0, ymax = panel_max)

figure_s4b <- ggplot(plot_data, aes(x = x, y = Pathogen_count, fill = Category)) +
  geom_rect(
    data = panel_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = Category),
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 0.45
  ) +
  geom_col(width = 0.73, color = NA, alpha = 0.92) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = category_fill, drop = FALSE) +
  scale_color_manual(values = category_border, guide = "none", drop = FALSE) +
  scale_x_continuous(
    limits = c(0.45, 8.55),
    breaks = seq_along(tissue_order),
    labels = rep("", length(tissue_order)),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Number of pathogens", tag = "b") +
  theme_bw(base_family = "Arial", base_size = 9.5) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#D9D9D9", linewidth = 0.35),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 8.5),
    axis.title.y = element_text(color = "black", size = 10.5, margin = margin(r = 5)),
    axis.ticks.y = element_line(color = "#777777", linewidth = 0.35),
    panel.spacing = grid::unit(0.12, "cm"),
    plot.tag = element_text(face = "bold", size = 17),
    plot.tag.position = c(0.005, 0.995),
    plot.margin = margin(t = 5, r = 5, b = 4, l = 5),
    legend.position = "none"
  )

ggsave(file.path(output_dir, "FigureS4B_pathogen_count_by_category.pdf"), figure_s4b, width = 4.2, height = 5.4, units = "in", device = grDevices::cairo_pdf, bg = "white")
ggsave(file.path(output_dir, "FigureS4B_pathogen_count_by_category.png"), figure_s4b, width = 4.2, height = 5.4, units = "in", dpi = 600, bg = "white")
