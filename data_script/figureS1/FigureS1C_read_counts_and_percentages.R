#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
required <- c("dplyr", "tidyr", "ggplot2", "patchwork", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

qc_filtered_colour <- "#BFC3C6"
rrna_filtered_colour <- "#AEA1E2"
host_mapped_colour <- "#C1B77E"
remaining_colour <- "#CFE2BF"
figure_width <- 14
figure_height <- 10.8
label_angle <- 55
label_size <- 2.5

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_paths <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) "" else as.character(frame$ofile)
}, character(1))
source_paths <- source_paths[nzchar(source_paths)]
editor_path <- ""
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  editor_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                          error = function(e) "")
}
script_dir <- if (length(script_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/"))
} else if (length(source_paths) > 0) {
  dirname(normalizePath(tail(source_paths, 1), winslash = "/"))
} else if (nzchar(editor_path)) {
  dirname(normalizePath(editor_path, winslash = "/"))
} else {
  getwd()
}
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else file.path(script_dir, "FigureS1C_library_read_qc.tsv")
output_prefix <- if (length(args) >= 2) args[2] else file.path(script_dir, "FigureS1C_read_counts_and_percentages")

qc <- utils::read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c(
  "library_id", "host_order", "order_position", "host_species", "total_reads",
  "after_qc_reads", "after_simple_filter_reads", "non_rrna_input_reads",
  "concordant_host_reads", "remaining_reads"
)
if (length(setdiff(required_columns, names(qc)))) stop("Incomplete Figure S1C input columns.")
if (nrow(qc) != 1027 || anyDuplicated(qc$library_id)) stop("Expected 1027 unique libraries.")
count_columns <- c("total_reads", "after_qc_reads", "after_simple_filter_reads",
                   "non_rrna_input_reads", "concordant_host_reads", "remaining_reads")
count_matrix <- as.matrix(qc[count_columns])
if (any(!is.finite(count_matrix)) || any(count_matrix < 0) || any(qc$total_reads <= 0)) {
  stop("Read counts must be finite and nonnegative, with positive total reads.")
}
if (any(qc$total_reads < qc$after_qc_reads) ||
    any(qc$after_qc_reads < qc$after_simple_filter_reads) ||
    any(qc$after_simple_filter_reads < qc$non_rrna_input_reads) ||
    any(qc$non_rrna_input_reads < qc$remaining_reads) ||
    any(qc$concordant_host_reads + qc$remaining_reads != qc$non_rrna_input_reads)) {
  stop("Read counts disagree with the retained processing stages.")
}

qc <- qc %>% mutate(
  qc_removed_reads = total_reads - after_simple_filter_reads,
  rrna_filtering_reads = after_simple_filter_reads - non_rrna_input_reads
)
orders <- qc %>% distinct(host_order, order_position) %>% arrange(order_position)
if (anyDuplicated(orders$host_order) || anyDuplicated(orders$order_position)) {
  stop("Host orders must have one unique plotting position each.")
}
qc <- qc %>%
  mutate(host_order = factor(host_order, levels = orders$host_order)) %>%
  arrange(host_order, desc(remaining_reads), library_id) %>%
  mutate(library_position = row_number())
order_centres <- qc %>% group_by(host_order) %>% summarise(
  n_libraries = n(), start = min(library_position), end = max(library_position),
  centre = (start + end) / 2, .groups = "drop"
)
boundaries <- order_centres %>% filter(end < nrow(qc)) %>% mutate(boundary = end + 0.5)

read_labels <- c(
  qc_removed_reads = "QC-filtered reads",
  rrna_filtering_reads = "Further filtering/rRNA reads",
  concordant_host_reads = "Host-mapped reads",
  remaining_reads = "Remaining reads"
)
read_colors <- setNames(c(qc_filtered_colour, rrna_filtered_colour,
                          host_mapped_colour, remaining_colour), unname(read_labels))
long_data <- qc %>%
  select(library_id, library_position, host_order, host_species, total_reads,
         all_of(names(read_labels))) %>%
  pivot_longer(cols = all_of(names(read_labels)), names_to = "read_class", values_to = "read_count") %>%
  mutate(read_class = factor(unname(read_labels[read_class]), levels = unname(read_labels)),
         percentage = 100 * read_count / total_reads,
         log10_total_reads = log10(total_reads),
         upper_segment_height = log10_total_reads * read_count / total_reads)
totals <- long_data %>% group_by(library_id) %>% summarise(
  count_sum = sum(read_count), expected = first(total_reads),
  percentage_sum = sum(percentage), .groups = "drop"
)
if (any(totals$count_sum != totals$expected) || any(abs(totals$percentage_sum - 100) > 1e-8)) {
  stop("Stacked counts or percentages failed the per-library sum check.")
}
utils::write.table(long_data, paste0(output_prefix, "_plot_data.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE, na = "")

if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))
common_theme <- theme_classic(base_family = "Arial", base_size = 10) + theme(
  axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  axis.text.y = element_text(size = 9), axis.title.y = element_text(size = 11),
  legend.title = element_text(face = "bold", size = 10),
  legend.text = element_text(size = 9),
  plot.tag = element_text(face = "bold", size = 16)
)

p_counts <- ggplot(long_data, aes(library_position, upper_segment_height, fill = read_class)) +
  geom_col(width = 1, colour = NA) +
  geom_vline(data = boundaries, aes(xintercept = boundary), colour = "white", linewidth = 0.6) +
  scale_fill_manual(values = read_colors, breaks = unname(read_labels), name = "Read type") +
  scale_x_continuous(limits = c(0.5, nrow(qc) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(
    breaks = sort(unique(c(scales::pretty_breaks(5)(c(0, max(long_data$log10_total_reads))), 9))),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(x = NULL, y = expression(log[10](total~reads)), tag = "c") + common_theme +
  theme(plot.margin = margin(8, 8, 8, 8))

p_percent <- ggplot(long_data, aes(library_position, percentage, fill = read_class)) +
  geom_col(width = 1, colour = NA) +
  geom_vline(data = boundaries, aes(xintercept = boundary), colour = "white", linewidth = 0.6) +
  scale_fill_manual(values = read_colors, breaks = unname(read_labels), name = "Read type") +
  scale_x_continuous(limits = c(0.5, nrow(qc) + 0.5),
                     breaks = order_centres$centre, labels = as.character(order_centres$host_order),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  labs(x = NULL, y = "Percentage of total reads") + common_theme +
  theme(
    axis.text.x = element_text(family = "Arial", colour = "black",
                               size = label_size * 72.27 / 25.4, angle = label_angle,
                               hjust = 1, vjust = 1),
    axis.ticks.x = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )

figure <- (p_counts / p_percent + plot_layout(heights = c(1, 1), guides = "collect")) &
  theme(legend.position = "top")
grDevices::pdfFonts(Arial = grDevices::pdfFonts("ArialMT")[[1]])
ggsave(paste0(output_prefix, ".pdf"), figure, width = figure_width, height = figure_height,
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(paste0(output_prefix, ".png"), figure, width = figure_width, height = figure_height, dpi = 600, bg = "white")
message("Saved Figure S1C for ", nrow(qc), " libraries in ", nrow(order_centres), " host orders.")
