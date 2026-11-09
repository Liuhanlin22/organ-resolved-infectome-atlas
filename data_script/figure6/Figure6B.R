#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
required <- c("openxlsx", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(ggplot2))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_paths <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) "" else as.character(frame$ofile)
}, character(1))
source_paths <- source_paths[nzchar(source_paths)]
script_dir <- if (length(script_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/"))
} else if (length(source_paths)) {
  dirname(normalizePath(tail(source_paths, 1), winslash = "/"))
} else {
  getwd()
}
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else file.path(script_dir, "Figure6B_input.xlsx")
output_prefix <- if (length(args) >= 2) args[2] else file.path(script_dir, "Figure6B_host_order_breadth")
if (!file.exists(input_file)) stop("Input not found: ", input_file)

dat <- openxlsx::read.xlsx(input_file, sheet = "Selected_GE2", check.names = FALSE, sep.names = " ")
required_columns <- c("Rank", "Updated pathogen name", "Host range (order)", "Host orders")
if (!all(required_columns %in% names(dat))) stop("Missing required Figure6B input columns.")
plot_data <- data.frame(
  rank = as.numeric(dat[["Rank"]]),
  pathogen = trimws(as.character(dat[["Updated pathogen name"]])),
  host_order_count = as.numeric(dat[["Host range (order)"]]),
  host_orders = trimws(as.character(dat[["Host orders"]]))
)
if (!nrow(plot_data) || anyNA(plot_data) || any(!nzchar(plot_data$pathogen)) ||
    anyDuplicated(plot_data$pathogen) || anyDuplicated(plot_data$rank)) {
  stop("Figure6B requires unique, nonempty taxa and ranks.")
}
if (any(!is.finite(plot_data$host_order_count)) ||
    any(plot_data$host_order_count < 2) ||
    any(plot_data$host_order_count != floor(plot_data$host_order_count))) {
  stop("Host-order counts must be integers of at least two.")
}
listed_counts <- vapply(strsplit(plot_data$host_orders, ";", fixed = TRUE), function(x) {
  length(unique(trimws(x[nzchar(trimws(x))])))
}, integer(1))
if (any(listed_counts != plot_data$host_order_count)) {
  stop("Host-order counts do not match the listed host orders.")
}
plot_data <- plot_data[order(plot_data$rank), ]
plot_data$pathogen <- factor(plot_data$pathogen, levels = plot_data$pathogen)
if (.Platform$OS.type == "windows") grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))

plot <- ggplot(plot_data, aes(x = pathogen, xend = pathogen)) +
  geom_segment(aes(y = 0, yend = host_order_count),
    colour = "#87E5EC", linewidth = 4.8, lineend = "round") +
  scale_x_discrete(position = "top", expand = expansion(add = 0.65)) +
  scale_y_continuous(breaks = seq(0, max(plot_data$host_order_count)),
    expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, max(plot_data$host_order_count) + 0.35)) +
  labs(x = NULL, y = "Host range (order)", tag = "b") +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5,
      size = 8, face = "italic", colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 11),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    plot.tag = element_text(size = 17, face = "bold"),
    plot.margin = margin(10, 8, 8, 8)
  )
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
ggsave(paste0(output_prefix, ".pdf"), plot, width = 5.5, height = 5.5,
  units = "in", device = pdf_device, bg = "white")
ggsave(paste0(output_prefix, ".png"), plot, width = 5.5, height = 5.5,
  units = "in", dpi = 600, bg = "white")
message("Figure6B complete: ", nrow(plot_data), " taxa, host-order counts ",
  min(plot_data$host_order_count), " to ", max(plot_data$host_order_count), ".")
