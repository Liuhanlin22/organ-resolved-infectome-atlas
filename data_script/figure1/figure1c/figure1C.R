#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
required <- c("dplyr", "ggplot2", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_paths <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) "" else as.character(frame$ofile)
}, character(1))
source_paths <- source_paths[nzchar(source_paths)]
script_dir <- if (length(script_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/"))
} else if (length(source_paths) > 0) {
  dirname(normalizePath(tail(source_paths, 1), winslash = "/"))
} else {
  getwd()
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else file.path(script_dir, "figure1C_species_library_count.tsv")
output_prefix <- if (length(args) >= 2) args[2] else file.path(script_dir, "figure1C_host_species_library_count")
dat <- utils::read.delim(
  input_file, check.names = FALSE, stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)
required_columns <- c("host_order", "species", "species_abbreviation", "library_count", "species_position")
if (length(setdiff(required_columns, names(dat)))) stop("Incomplete Figure 1C input columns.")
dat <- dplyr::arrange(dat, species_position)
if (nrow(dat) != 51 || sum(dat$library_count) != 1027 || anyDuplicated(dat$species)) {
  stop("Expected 51 unique species and 1027 libraries.")
}
order_levels <- unique(dat$host_order)
dat <- dat %>% mutate(
  host_order = factor(host_order, levels = order_levels),
  species = factor(species, levels = species)
)
order_palette_all <- c(
  "Anguilliformes" = "#54278F", "Aulopiformes" = "#756BB1",
  "Carangaria incertae sedis" = "#9E9AC8", "Carangiformes" = "#BCBDDC",
  "Chaetodontiformes" = "#DADAEB", "Clupeiformes" = "#1C69AF",
  "Eupercaria incertae sedis" = "#3182BD", "Gadiformes" = "#6BAED6",
  "Gobiiformes" = "#9ECAE1", "Lutjaniformes" = "#C6DBEF",
  "Mugiliformes" = "#006D2C", "Myctophiformes" = "#31A354",
  "Perciformes" = "#74C476", "Pleuronectiformes" = "#A1D99B",
  "Scombriformes" = "#C7E9C0", "Siluriformes" = "#C60505",
  "Spariformes" = "#DE2D26", "Tetraodontiformes" = "#FC9272",
  "Carcharhiniformes" = "#FCBBA1", "Myliobatiformes" = "#F4D4C9"
)
if (length(setdiff(order_levels, names(order_palette_all)))) stop("Missing host-order palette entry.")
boundaries <- dat %>% count(host_order, name = "n_species") %>%
  mutate(boundary = cumsum(n_species) + 0.5) %>% filter(row_number() < n())
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))
plot <- ggplot(dat, aes(x = species, y = library_count, fill = host_order)) +
  geom_col(width = 0.82) +
  geom_vline(data = boundaries, aes(xintercept = boundary), linetype = "dashed",
             linewidth = 0.25, colour = "#999999") +
  scale_fill_manual(values = order_palette_all[order_levels], name = "Host order/order-level lineage") +
  scale_x_discrete(labels = setNames(dat$species_abbreviation, as.character(dat$species))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06)), breaks = scales::pretty_breaks(5)) +
  labs(x = NULL, y = "Host species library count", tag = "c") +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_classic(base_family = "Arial") +
  theme(
    legend.position = "top", legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 6.5),
    axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1, size = 6.5, face = "italic"),
    axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 10),
    plot.tag = element_text(face = "bold", size = 16), plot.margin = margin(8, 8, 8, 8)
  )
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
ggsave(paste0(output_prefix, ".pdf"), plot, width = 13.5, height = 5.5, device = pdf_device)
ggsave(paste0(output_prefix, ".png"), plot, width = 13.5, height = 5.5, dpi = 600, bg = "white")
