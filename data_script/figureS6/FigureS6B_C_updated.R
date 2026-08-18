#!/usr/bin/env Rscript
required <- c("openxlsx", "ggplot2", "grid")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(openxlsx); library(ggplot2); library(grid)})

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(arg) == 1) dirname(normalizePath(sub("^--file=", "", arg), winslash = "/")) else normalizePath(getwd(), winslash = "/")
input_file <- file.path(script_dir, "FigureS6B_C_input_updated.xlsx")
all_data <- read.xlsx(input_file, sheet = "all_75_range", check.names = FALSE)

stopifnot(nrow(all_data) == 75)
stopifnot(!anyDuplicated(all_data$pathogen))
if (any(!is.finite(all_data$total_RPM)) || any(all_data$total_RPM < 0)) stop("Invalid total_RPM values")

range_order <- c("Within species", "Cross species", "Cross genus", "Cross family", "Cross order")
category_order <- c("RNA", "DNA", "BACTERIA", "EUKARYOTA")
taxa_colours <- c(
  "Mono-Chu" = "#b9d7f4", "Reo" = "#f7aaa6", "Noda" = "#b8a8dc",
  "Bacillota" = "#9ed3e6", "Pseudomonadota" = "#39c1c7",
  "Nematoda" = "#ee5b8f", "Platyhelminthes" = "#b00076"
)
category_colours <- c(
  RNA = "#F2DCDB", DNA = "#B8CAE0", BACTERIA = "#C8B7D9", EUKARYOTA = "#CFE0B6"
)

b_data <- all_data[all_data$cross_range_code > 0, , drop = FALSE]
b_data$cross_range <- factor(b_data$cross_range, levels = range_order[2:5])
b_data$taxa_display <- factor(b_data$taxa_display, levels = names(taxa_colours))
b_data$pathogen <- factor(b_data$pathogen, levels = rev(as.character(b_data$pathogen)))

counts <- as.data.frame(table(
  factor(all_data$cross_range, levels = range_order),
  factor(all_data$category_4way, levels = category_order),
  dnn = c("cross_range", "category")
))
names(counts)[3] <- "n"
counts$category <- factor(counts$category, levels = category_order)
counts$cross_range <- factor(counts$cross_range, levels = range_order)
if (sum(counts$n) != 75) stop("FigureS6C total is not 75")

p_b <- ggplot(b_data, aes(x = cross_range, y = pathogen, colour = taxa_display)) +
  geom_point(size = 4.4, show.legend = TRUE) +
  scale_x_discrete(drop = FALSE, limits = range_order[2:5]) +
  scale_colour_manual(values = taxa_colours, drop = TRUE, name = "Taxa") +
  labs(tag = "b", x = NULL, y = NULL) +
  theme_classic(base_size = 9, base_family = "sans") +
  theme(
    plot.tag = element_text(face = "bold", size = 14, hjust = 0, vjust = 1),
    axis.text.y = element_text(face = "italic", size = 8),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 8),
    axis.ticks = element_line(linewidth = 0.35),
    panel.grid.major.y = element_line(colour = "#e8e8e8", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key = element_rect(fill = NA, colour = NA),
    plot.margin = margin(4, 4, 4, 4)
  )

p_c <- ggplot(counts, aes(x = cross_range, y = n, fill = category)) +
  geom_col(width = 0.58, colour = NA, position = position_stack(reverse = TRUE)) +
  scale_x_discrete(drop = FALSE, limits = range_order) +
  scale_fill_manual(values = category_colours, breaks = category_order, drop = FALSE, name = "Category") +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 10), expand = c(0, 0)) +
  labs(tag = "c", x = NULL, y = "Number of pathogens") +
  theme_classic(base_size = 9, base_family = "sans") +
  theme(
    plot.tag = element_text(face = "bold", size = 14, hjust = 0, vjust = 1),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 8),
    axis.ticks = element_line(linewidth = 0.35),
    panel.grid.major.y = element_line(colour = "#d9d9d9", linetype = "dashed", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key = element_rect(fill = NA, colour = NA),
    plot.margin = margin(4, 4, 4, 4)
  )

save_plot <- function(plot, base, width, height) {
  ggsave(paste0(base, ".png"), plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
  ggsave(paste0(base, ".tiff"), plot, width = width, height = height, units = "in", dpi = 300, compression = "lzw", bg = "white")
  grDevices::pdf(paste0(base, ".pdf"), width = width, height = height,
                 useDingbats = FALSE, compress = TRUE)
  print(plot)
  dev.off()
  svg(paste0(base, ".svg"), width = width, height = height, bg = "white")
  print(plot)
  dev.off()
}

save_plot(p_b, file.path(script_dir, "FigureS6B_R_rendered"), 5.6, 4.6)
save_plot(p_c, file.path(script_dir, "FigureS6C_R_rendered"), 4.8, 4.2)

# Combined panel export without requiring patchwork or gridExtra.
pdf(file.path(script_dir, "FigureS6B_C_R_rendered.pdf"), width = 10.8, height = 4.6, useDingbats = FALSE)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1.35, 1), "null"))))
print(p_b, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_c, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
popViewport()
dev.off()

png(file.path(script_dir, "FigureS6B_C_R_rendered.png"), width = 3240, height = 1380, res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1.35, 1), "null"))))
print(p_b, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_c, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
popViewport()
dev.off()
