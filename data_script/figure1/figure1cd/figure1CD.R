#!/usr/bin/env Rscript

# Figure 1C-D — host library coverage and read-quality summary
# Input: Supplement_Table1_Info_sample, related to figure1-revised.xlsx
# Run: Rscript figure1CD.R "input.xlsx" "figure1CD"
# Outputs: species-library table, library-QC table, PDF and PNG figure.

options(stringsAsFactors = FALSE)
required_packages <- c("readxl", "dplyr", "tidyr", "ggplot2", "patchwork", "scales", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) install.packages(missing_packages, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(scales); library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "Supplement_Table1_Info_sample, related to figure1-revised.xlsx"
output_prefix <- if (length(args) >= 2) args[2] else "figure1CD"
if (!file.exists(input_file)) stop("Input Excel file not found: ", input_file)

# Appearance parameters -------------------------------------------------------
font_family <- "Arial"
species_label_size <- 6.5
order_label_size <- 6.5
fig_width <- 13.5
fig_height <- 11.0
rRNA_colour <- "#92A9D0"
remaining_colour <- "#DDEBD3"

# Register Arial explicitly on Windows, then use Cairo for PDF export so the
# requested system font is handled consistently in both PDF and PNG outputs.
if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

# Locate the true header row. The supplied workbook has a descriptive title in
# row 1 and the actual column names in row 2; automatic detection is safer.
sheet_name <- excel_sheets(input_file)[1]
preview <- read_excel(input_file, sheet = sheet_name, col_names = FALSE, n_max = 10, .name_repair = "minimal")
header_row <- which(apply(preview, 1, function(z) any(trimws(as.character(z)) == "lib_index", na.rm = TRUE)))[1]
if (is.na(header_row)) stop("Could not locate the header row containing 'lib_index'.")
raw <- read_excel(input_file, sheet = sheet_name, skip = header_row - 1, .name_repair = "minimal")
names(raw) <- trimws(names(raw))

required_columns <- c("lib_index", "order / order-level lineage", "species", "species abbreviation", "totol_reads", "after_qc_reads", "remove_complicated_region_reads", "remove_rRNA_reads", "remaining reads (%)")
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns) > 0) stop("Required input column(s) missing: ", paste(missing_columns, collapse = ", "))

# Figure order: bony-fish orders alphabetically from left to right, followed by
# the two cartilaginous-fish orders at the far right.
order_levels_requested <- c(
  "Anguilliformes", "Aulopiformes", "Carangaria incertae sedis",
  "Carangiformes", "Chaetodontiformes", "Clupeiformes",
  "Eupercaria incertae sedis", "Gadiformes", "Gobiiformes",
  "Lutjaniformes", "Mugiliformes", "Myctophiformes", "Perciformes",
  "Pleuronectiformes", "Scombriformes", "Siluriformes", "Spariformes",
  "Tetraodontiformes", "Carcharhiniformes", "Myliobatiformes"
)
orders_observed <- raw %>% transmute(host_order = trimws(`order / order-level lineage`)) %>% filter(!is.na(host_order), host_order != "") %>% distinct() %>% pull(host_order)
order_levels <- c(order_levels_requested[order_levels_requested %in% orders_observed],
                  sort(setdiff(orders_observed, order_levels_requested)))
dat <- raw %>% transmute(
  lib_index = trimws(as.character(lib_index)),
  host_order = trimws(`order / order-level lineage`),
  species = trimws(species), species_abbreviation = trimws(`species abbreviation`),
  total_reads = as.numeric(totol_reads), qc_filtered_reads = as.numeric(after_qc_reads),
  pre_rRNA_reads = as.numeric(remove_complicated_region_reads),
  remaining_reads = as.numeric(remove_rRNA_reads),
  remaining_reads_pct_input = as.numeric(`remaining reads (%)`)
) %>% filter(!is.na(lib_index), lib_index != "", !is.na(host_order), host_order != "", !is.na(species), species != "") %>% mutate(
  host_order = factor(host_order, levels = order_levels),
  species_abbreviation = if_else(is.na(species_abbreviation) | species_abbreviation == "", species, species_abbreviation),
  rRNA_reads = pmax(pre_rRNA_reads - remaining_reads, 0),
  # All reads removed before the rRNA step: QC loss plus simple/complicated
  # region filtering.  Together with rRNA_reads and remaining_reads these sum
  # exactly to total_reads, allowing a true 100%-stacked three-colour bar.
  qc_filtered_loss_reads = pmax(total_reads - pre_rRNA_reads, 0),
  remaining_reads_pct = if_else(pre_rRNA_reads > 0, 100 * remaining_reads / pre_rRNA_reads, NA_real_),
  rRNA_reads_pct = if_else(pre_rRNA_reads > 0, 100 * rRNA_reads / pre_rRNA_reads, NA_real_),
  qc_retained_pct = if_else(total_reads > 0, 100 * qc_filtered_reads / total_reads, NA_real_),
  qc_filtered_loss_pct = if_else(total_reads > 0, 100 * qc_filtered_loss_reads / total_reads, NA_real_),
  rRNA_reads_pct_total = if_else(total_reads > 0, 100 * rRNA_reads / total_reads, NA_real_),
  remaining_reads_pct_total = if_else(total_reads > 0, 100 * remaining_reads / total_reads, NA_real_)
)
if (any(is.na(dat$pre_rRNA_reads)) || any(is.na(dat$remaining_reads))) stop("Missing values in rRNA filtering columns.")
if (any(abs(dat$remaining_reads_pct + dat$rRNA_reads_pct - 100) > 0.001, na.rm = TRUE)) stop("rRNA and remaining proportions do not sum to 100%.")

# Two processed tables requested for Figure 1C and Figure 1D.
species_library_count <- dat %>% group_by(host_order, species, species_abbreviation) %>% summarise(library_count = n(), .groups = "drop") %>% arrange(host_order, species_abbreviation)
library_read_qc <- dat %>% select(host_order, species, species_abbreviation, lib_index, total_reads, qc_filtered_reads, qc_retained_pct, pre_rRNA_reads, qc_filtered_loss_reads, qc_filtered_loss_pct, rRNA_reads, rRNA_reads_pct, rRNA_reads_pct_total, remaining_reads, remaining_reads_pct, remaining_reads_pct_total, remaining_reads_pct_input) %>% arrange(host_order, species_abbreviation, lib_index)
write_tsv(species_library_count, paste0(output_prefix, "_species_library_count.tsv"))
write_tsv(library_read_qc, paste0(output_prefix, "_library_read_qc.tsv"))

order_levels_used <- levels(droplevels(dat$host_order))
# Fixed palette matching the previously used Figure 1C order colours.  The
# entries not present in the older panel (Eupercaria, Chaetodontiformes,
# Carangaria incertae sedis and Lutjaniformes) are explicitly assigned here.
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
unknown_orders <- setdiff(order_levels_used, names(order_palette_all))
if (length(unknown_orders) > 0) {
  order_palette_all[unknown_orders] <- grDevices::hcl(
    h = seq(15, 345, length.out = length(unknown_orders)), c = 60, l = 58
  )
  warning("New host order(s) received automatically assigned colours: ", paste(unknown_orders, collapse = ", "))
}
order_palette <- order_palette_all[order_levels_used]
species_library_count <- species_library_count %>% mutate(
  species_abbreviation = factor(species_abbreviation, levels = unique(species_abbreviation)),
  host_order = factor(host_order, levels = order_levels_used)
)

# Figure 1C: number of libraries available for each host species.
c_boundaries <- species_library_count %>% count(host_order, name = "n_species") %>%
  mutate(boundary = cumsum(n_species) + 0.5) %>% filter(row_number() < n())
p_c <- ggplot(species_library_count, aes(x = species_abbreviation, y = library_count, fill = host_order)) +
  geom_col(width = 0.82) +
  geom_vline(data = c_boundaries, aes(xintercept = boundary), linetype = "dashed", linewidth = 0.25, colour = "#999999") +
  scale_fill_manual(values = order_palette, name = "Host order") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06)), breaks = pretty_breaks(5)) +
  labs(x = NULL, y = "Host species library count", tag = "c") + theme_classic(base_family = font_family) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "top", legend.title = element_text(face = "bold", size = 9), legend.text = element_text(size = order_label_size), axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1, size = species_label_size), axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 10), plot.tag = element_text(face = "bold", size = 16), plot.tag.position = c(-0.03, 1.03), plot.margin = margin(5.5, 5.5, 2, 5.5))

# Figure 1D: three mutually exclusive read categories for every library.
library_read_qc <- library_read_qc %>% mutate(host_order = factor(host_order, levels = order_levels_used), species_abbreviation = factor(species_abbreviation, levels = levels(species_library_count$species_abbreviation))) %>% arrange(host_order, species_abbreviation, lib_index) %>% mutate(library_num = row_number())
d_boundaries <- library_read_qc %>% count(host_order, name = "n_libraries") %>%
  mutate(boundary = cumsum(n_libraries) + 0.5) %>% filter(row_number() < n())
order_centres <- library_read_qc %>% count(host_order, name = "n_libraries") %>%
  mutate(end = cumsum(n_libraries), start = lag(end, default = 0) + 1,
         centre = (start + end) / 2,
         # Wide order groups can be labelled near 0%; narrow groups receive a
         # longer individual connector and are placed farther below.
         label_y = case_when(
           n_libraries >= 80 ~ -7,
           n_libraries >= 40 ~ -14,
           n_libraries >= 20 ~ -22,
           TRUE ~ -28 - 5 * ((row_number() - 1) %% 3)
         ),
         connector_end = label_y + 2.5)
d_long <- library_read_qc %>% select(library_num, qc_filtered_loss_pct, rRNA_reads_pct_total, remaining_reads_pct_total) %>% pivot_longer(cols = c(qc_filtered_loss_pct, rRNA_reads_pct_total, remaining_reads_pct_total), names_to = "read_class", values_to = "percentage") %>% mutate(read_class = recode(read_class, qc_filtered_loss_pct = "QC-filtered reads", rRNA_reads_pct_total = "rRNA reads", remaining_reads_pct_total = "Remaining reads"), read_class = factor(read_class, levels = c("QC-filtered reads", "rRNA reads", "Remaining reads"))) %>% group_by(library_num) %>% mutate(percentage = 100 * percentage / sum(percentage)) %>% ungroup()
p_d <- ggplot(d_long, aes(x = library_num, y = percentage, fill = read_class)) +
  geom_col(width = 1, colour = NA) +
  # White gaps, rather than dark internal lines, distinguish order groups.
  geom_vline(data = d_boundaries, aes(xintercept = boundary), colour = "white", linewidth = 1.1) +
  # One connector per order, directly centred on its own library block.
  geom_segment(data = order_centres, aes(x = centre, xend = centre, y = 0, yend = connector_end), inherit.aes = FALSE, colour = "#444444", linewidth = 0.35) +
  geom_text(data = order_centres, aes(x = centre, y = label_y, label = host_order), inherit.aes = FALSE, size = order_label_size / 2.8, family = font_family, angle = 45, hjust = 1, vjust = 1) +
  # Draw the y axis explicitly so it ends exactly at 0%, not beneath the bars.
  annotate("segment", x = 0.5, xend = 0.5, y = 0, yend = 100, linewidth = 0.55, colour = "black") +
  scale_fill_manual(values = c("QC-filtered reads" = "#E7B66B", "rRNA reads" = rRNA_colour, "Remaining reads" = remaining_colour)) +
  # Lock the first and last library to the panel borders.  This removes the
  # unwanted left/right blank space and matches Figure 1C's plot width.
  scale_x_continuous(limits = c(0.5, nrow(library_read_qc) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  coord_cartesian(ylim = c(-45, 100), clip = "off") +
  labs(x = NULL, y = "Reads proportion", fill = NULL, tag = "d") + theme_classic(base_family = font_family) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 10), legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 9), plot.tag = element_text(face = "bold", size = 16), plot.tag.position = c(-0.03, 1.03), plot.margin = margin(5.5, 5.5, 47, 5.5))

final_figure <- p_c / p_d + plot_layout(heights = c(1.05, 1))
# cairo_pdf accesses the Windows system-font registry; it avoids repeated Arial
# warnings produced by the basic PostScript PDF device.
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
ggsave(paste0(output_prefix, ".pdf"), final_figure, width = fig_width, height = fig_height, units = "in", device = pdf_device)
ggsave(paste0(output_prefix, ".png"), final_figure, width = fig_width, height = fig_height, units = "in", dpi = 600, bg = "white")
