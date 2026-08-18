#!/usr/bin/env Rscript

# Figure 5G — pathogen association coefficients across 112 fish
# Statistics: phi coefficient for binary detection profiles; one-sided
# Fisher exact tests; Benjamini-Hochberg correction across all detected pairs.

required_packages <- c("openxlsx", "writexl", "ggplot2", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}

input_file <- file.path(script_dir, "Figure5G_input.xlsx")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

rpm_df <- read.xlsx(input_file, sheet = "host_by_pathogen_RPM", check.names = FALSE)
annotations <- read.xlsx(input_file, sheet = "pathogen_annotations", check.names = FALSE)

if (!"host_id" %in% names(rpm_df)) stop("Missing host_id column")
matrix_columns <- setdiff(names(rpm_df), "host_id")
pathogen_names <- annotations$Species.name
if (nrow(rpm_df) != 112 || length(matrix_columns) != 75 || length(pathogen_names) != 75) {
  stop("Expected 112 fish x 75 pathogens")
}
if (anyDuplicated(rpm_df$host_id) || anyDuplicated(pathogen_names)) {
  stop("Duplicated host or pathogen identifiers detected")
}

rpm <- as.matrix(rpm_df[, matrix_columns, drop = FALSE])
storage.mode(rpm) <- "numeric"
colnames(rpm) <- pathogen_names
if (anyNA(rpm) || any(rpm < 0)) stop("RPM matrix contains missing or negative values")

present <- rpm > 0
infected_hosts <- colSums(present)
detected <- pathogen_names[infected_hosts > 0]
displayed <- pathogen_names[infected_hosts >= 3]

if (length(detected) != 68 || length(displayed) != 16) {
  stop(
    "Expected 68 detected and 16 displayed pathogens; observed ",
    length(detected), " and ", length(displayed)
  )
}

phi_matrix <- cor(present[, detected, drop = FALSE] * 1, method = "pearson")
diag(phi_matrix) <- 1

pair_index <- combn(detected, 2, simplify = FALSE)
pairwise <- do.call(
  rbind,
  lapply(pair_index, function(pair) {
    x <- present[, pair[1]]
    y <- present[, pair[2]]
    both <- sum(x & y)
    x_only <- sum(x & !y)
    y_only <- sum(!x & y)
    neither <- sum(!x & !y)
    ft <- fisher.test(
      matrix(c(both, x_only, y_only, neither), nrow = 2, byrow = TRUE),
      alternative = "greater"
    )
    data.frame(
      pathogen_1 = pair[1],
      pathogen_2 = pair[2],
      infected_hosts_1 = sum(x),
      infected_hosts_2 = sum(y),
      co_infected_hosts = both,
      expected_co_infected_hosts = sum(x) * sum(y) / length(x),
      phi_coefficient = phi_matrix[pair[1], pair[2]],
      odds_ratio = unname(ft$estimate),
      fisher_p = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
)
pairwise$fisher_BH_q <- p.adjust(pairwise$fisher_p, method = "BH")
pairwise$significant_positive_BH <- with(
  pairwise,
  fisher_BH_q < 0.05 & phi_coefficient > 0
)
pairwise <- pairwise[order(pairwise$fisher_BH_q, pairwise$fisher_p), ]

display_phi <- phi_matrix[displayed, displayed, drop = FALSE]
hc <- hclust(as.dist(1 - display_phi), method = "complete")
display_order <- rownames(display_phi)[hc$order]

display_info <- data.frame(
  display_order = seq_along(display_order),
  pathogen = display_order,
  infected_hosts = as.integer(infected_hosts[display_order]),
  stringsAsFactors = FALSE
)

grid <- expand.grid(
  row_pathogen = display_order,
  col_pathogen = display_order,
  stringsAsFactors = FALSE
)
grid$row_index <- match(grid$row_pathogen, display_order)
grid$col_index <- match(grid$col_pathogen, display_order)
grid <- grid[grid$row_index > grid$col_index, ]

significant <- pairwise[
  pairwise$significant_positive_BH &
    pairwise$pathogen_1 %in% displayed &
    pairwise$pathogen_2 %in% displayed,
]

if (nrow(significant) > 0) {
  significant$index_1 <- match(significant$pathogen_1, display_order)
  significant$index_2 <- match(significant$pathogen_2, display_order)
  significant$col_pathogen <- ifelse(
    significant$index_1 < significant$index_2,
    significant$pathogen_1,
    significant$pathogen_2
  )
  significant$row_pathogen <- ifelse(
    significant$index_1 > significant$index_2,
    significant$pathogen_1,
    significant$pathogen_2
  )
  significant$label <- sprintf("%.2f", significant$phi_coefficient)
}

grid$col_pathogen <- factor(grid$col_pathogen, levels = display_order)
grid$row_pathogen <- factor(grid$row_pathogen, levels = rev(display_order))
if (nrow(significant) > 0) {
  significant$col_pathogen <- factor(significant$col_pathogen, levels = display_order)
  significant$row_pathogen <- factor(significant$row_pathogen, levels = rev(display_order))
}

p <- ggplot(grid, aes(x = col_pathogen, y = row_pathogen)) +
  geom_tile(fill = "white", colour = "#C7C7C7", linewidth = 0.32) +
  geom_point(
    data = significant,
    aes(size = phi_coefficient, fill = phi_coefficient),
    shape = 21, colour = "#9E2F52", stroke = 0.30
  ) +
  geom_text(
    data = significant,
    aes(label = label),
    family = "Arial", fontface = "bold", size = 3.1, colour = "black"
  ) +
  scale_size_continuous(range = c(5.5, 10.5), limits = c(0, 1), guide = "none") +
  scale_fill_gradient(
    name = "Correlation\ncoefficients",
    low = "#FFFFFF", high = "#B83B5E", limits = c(0, 1),
    breaks = seq(0, 1, 0.2), oob = squish
  ) +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  coord_fixed(clip = "off") +
  labs(tag = "g", x = NULL, y = NULL) +
  theme_void(base_family = "Arial", base_size = 10) +
  theme(
    axis.text.x.top = element_text(
      family = "Arial", face = "italic", colour = "black",
      angle = 90, hjust = 0, vjust = 0.5, size = 8.7,
      margin = margin(b = 4)
    ),
    axis.text.y = element_text(
      family = "Arial", face = "italic", colour = "black",
      size = 8.7, margin = margin(r = 5)
    ),
    legend.position = "right",
    legend.title = element_text(family = "Arial", size = 10),
    legend.text = element_text(family = "Arial", size = 8.5),
    plot.tag = element_text(family = "Arial", face = "bold", size = 20),
    plot.tag.position = c(-0.04, 1.08),
    plot.margin = margin(52, 12, 8, 14)
  ) +
  guides(fill = guide_colourbar(barheight = unit(35, "mm"), barwidth = unit(5, "mm")))

ggsave(
  file.path(script_dir, "Figure5G.pdf"), p,
  width = 7.4, height = 6.6, units = "in",
  device = cairo_pdf, bg = "white"
)
ggsave(
  file.path(script_dir, "Figure5G.png"), p,
  width = 7.4, height = 6.6, units = "in",
  dpi = 600, bg = "white"
)

write.csv(pairwise, file.path(script_dir, "Figure5G_pairwise_Fisher_BH.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(display_info, file.path(script_dir, "Figure5G_display_pathogens.csv"), row.names = FALSE, fileEncoding = "UTF-8")

summary_table <- data.frame(
  metric = c(
    "Total fish", "Positive fish", "Negative fish", "Total pathogens",
    "Detected pathogens tested", "Pairwise tests", "Displayed pathogens (>=3 hosts)",
    "BH-significant positive pairs", "BH-significant pairs visible in panel"
  ),
  value = c(
    nrow(present), sum(rowSums(present) > 0), sum(rowSums(present) == 0),
    ncol(present), length(detected), nrow(pairwise), length(displayed),
    sum(pairwise$significant_positive_BH), nrow(significant)
  )
)

phi_out <- data.frame(pathogen = rownames(phi_matrix), phi_matrix, check.names = FALSE)
display_phi_out <- data.frame(pathogen = rownames(display_phi), display_phi, check.names = FALSE)

writexl::write_xlsx(
  list(
    summary = summary_table,
    display_pathogens = display_info,
    significant_positive = significant,
    pairwise_Fisher_BH = pairwise,
    phi_matrix_68 = phi_out,
    phi_matrix_display16 = display_phi_out
  ),
  path = file.path(script_dir, "Figure5G_results.xlsx")
)
