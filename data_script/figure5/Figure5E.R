#!/usr/bin/env Rscript

# Figure 5E: individual-level pathogen-composition t-SNE
# Each row is one fish after summing its eight organ libraries.
# Only pathogen-positive fish are included (72 fish x 75 pathogens).

required <- c("readr", "dplyr", "ggplot2", "Rtsne", "vegan")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R package(s): ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(Rtsne)
  library(vegan)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}

rpm_file <- file.path(script_dir, "rpm_table_virus_masked.csv")
meta_file <- file.path(script_dir, "lib_meta.csv")
pathogen_file <- file.path(script_dir, "virus_meta_merged.csv")
stopifnot(file.exists(rpm_file), file.exists(meta_file), file.exists(pathogen_file))

rpm <- read_csv(rpm_file, show_col_types = FALSE)
meta <- read_csv(meta_file, show_col_types = FALSE)
pathogen_meta <- read_csv(pathogen_file, show_col_types = FALSE)

pathogens <- pathogen_meta$virus_name[as.logical(pathogen_meta$is_core)]
stopifnot(length(pathogens) == 75L, length(unique(pathogens)) == 75L)
stopifnot(all(c("lib_id", "order", "plot_group") %in% names(meta)))
stopifnot(all(c("lib_id", pathogens) %in% names(rpm)))

rpm <- rpm %>% select(lib_id, all_of(pathogens)) %>% arrange(lib_id)
meta <- meta %>% arrange(lib_id)
stopifnot(identical(rpm$lib_id, meta$lib_id))

x <- as.matrix(rpm[, pathogens])
storage.mode(x) <- "double"
stopifnot(nrow(x) == 72L, ncol(x) == 75L)
stopifnot(!anyNA(x), all(x >= 0), all(rowSums(x) > 0))
stopifnot(!any(duplicated(as.data.frame(x))))

# Fixed display scheme requested by the author. Current-order values outside
# these six orders are grouped as Others; taxonomy in lib_meta.csv is unchanged.
group_levels <- c(
  "Anguilliformes", "Aulopiformes", "Clupeiformes",
  "Gadiformes", "Gobiiformes", "Others", "Scorpaeniformes"
)
meta$plot_group <- ifelse(
  meta$order %in% setdiff(group_levels, "Others"), meta$order, "Others"
)
meta$plot_group <- factor(meta$plot_group, levels = group_levels)

group_colors <- c(
  "Anguilliformes" = "#67a8cd",
  "Aulopiformes" = "#414985",
  "Clupeiformes" = "#ff9d9f",
  "Gadiformes" = "#d43f3b",
  "Gobiiformes" = "#eea235",
  "Others" = "#9ed17b",
  "Scorpaeniformes" = "#832440"
)
present_groups <- group_levels[group_levels %in% as.character(unique(meta$plot_group))]

x_log <- log1p(x)
set.seed(20260814)
tsne <- Rtsne(
  x_log,
  dims = 2,
  perplexity = 6,
  max_iter = 5000,
  pca = TRUE,
  check_duplicates = TRUE,
  verbose = FALSE
)

plot_data <- meta %>%
  mutate(tSNE1 = tsne$Y[, 1], tSNE2 = tsne$Y[, 2])

# PERMANOVA is calculated from Bray-Curtis distances in the original
# log1p-abundance space, never from the t-SNE coordinates.
bray <- vegdist(x_log, method = "bray")
stats_group <- droplevels(meta$plot_group)
set.seed(20260814)
permanova <- adonis2(bray ~ stats_group, permutations = 9999)

# Test multivariate dispersion because unequal dispersion can affect PERMANOVA.
dispersion <- betadisper(
  bray,
  stats_group,
  type = "median",
  bias.adjust = TRUE,
  add = "lingoes"
)
set.seed(20260814)
permdisp <- permutest(dispersion, permutations = 9999)

r2 <- permanova$R2[1]
p_perm <- permanova$`Pr(>F)`[1]
p_disp <- permdisp$tab[1, "Pr(>F)"]
format_p <- function(p) if (p < 0.001) "< 0.001" else paste0("= ", formatC(p, digits = 3, format = "f"))

annotation <- paste0(
  "Bray-Curtis PERMANOVA\n",
  "R\u00b2 = ", sprintf("%.3f", r2), "\n",
  "p ", format_p(p_perm)
)

p <- ggplot(plot_data, aes(tSNE1, tSNE2, colour = plot_group)) +
  geom_point(shape = 16, size = 2.8) +
  scale_colour_manual(values = group_colors, breaks = present_groups, drop = TRUE) +
  annotate(
    "text",
    x = min(plot_data$tSNE1) + 0.025 * diff(range(plot_data$tSNE1)),
    y = min(plot_data$tSNE2) + 0.025 * diff(range(plot_data$tSNE2)),
    label = annotation,
    hjust = 0, vjust = 0, family = "Arial", size = 3.2
  ) +
  labs(
    tag = "e", x = "T-SNE axis 1", y = "T-SNE axis 2",
    colour = "Host order"
  ) +
  coord_equal(clip = "off") +
  theme_bw(base_family = "Arial", base_size = 9) +
  theme(
    panel.grid.major = element_line(colour = "#D9D9D9", linewidth = 0.35),
    panel.grid.minor = element_line(colour = "#EEEEEE", linewidth = 0.25),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.key = element_blank(),
    axis.text = element_text(colour = "black"),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 1),
    plot.margin = margin(7, 7, 7, 9)
  ) +
  guides(colour = guide_legend(title.position = "top", nrow = 2, byrow = TRUE))

group_counts <- plot_data %>% count(order, plot_group, name = "n_positive_fish")
write_csv(group_counts, file.path(script_dir, "Figure5E_group_counts.csv"))

statistics <- tibble(
  test = c("PERMANOVA", "PERMDISP"),
  statistic = c(permanova$F[1], permdisp$tab[1, "F"]),
  R2 = c(r2, NA_real_),
  p_value = c(p_perm, p_disp),
  permutations = 9999L
)
write_csv(statistics, file.path(script_dir, "Figure5E_statistics.csv"))

ggsave(
  file.path(script_dir, "Figure5E.pdf"), p,
  width = 6.6, height = 5.2, units = "in", device = grDevices::cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(script_dir, "Figure5E.png"), p,
  width = 6.6, height = 5.2, units = "in", dpi = 600, bg = "white"
)
