library(openxlsx)
library(pheatmap)
library(grid)

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()

input_xlsx <- file.path(script_dir, "viral_data_pheatmap.xlsx")
data_raw <- openxlsx::read.xlsx(input_xlsx, sheet = 1, check.names = FALSE)
if (!"lib_id" %in% colnames(data_raw)) stop("viral_data_pheatmap.xlsx must contain a lib_id column")

desired_order <- scan(file.path(script_dir, "order_data.txt"), what = "", sep = "\t", quiet = TRUE)
desired_order_row <- scan(file.path(script_dir, "order_data_row.txt"), what = "", sep = "\t", quiet = TRUE)
missing_cols <- setdiff(desired_order, colnames(data_raw))
if (length(missing_cols)) stop("Missing host columns: ", paste(missing_cols, collapse = ", "))
missing_rows <- setdiff(desired_order_row, data_raw$lib_id)
if (length(missing_rows)) stop("Missing pathogen rows: ", paste(missing_rows, collapse = ", "))

# Keep lib_id separate; the old script accidentally dropped the first host column.
row_ids <- as.character(data_raw$lib_id)
data <- data_raw[, desired_order, drop = FALSE]
rownames(data) <- row_ids
data <- data[desired_order_row, , drop = FALSE]
heatmap_data <- log10(as.matrix(data) + 1)

annotation_col <- read.delim(file.path(script_dir, "col_heatmap_annot.txt"),
                             header = TRUE, sep = "\t", row.names = 1,
                             check.names = FALSE, stringsAsFactors = FALSE,
                             quote = "")
annotation_row <- read.delim(file.path(script_dir, "row_heatmap_annot.txt"),
                             header = TRUE, sep = "\t", row.names = 1,
                             check.names = FALSE, stringsAsFactors = FALSE,
                             quote = "")
annotation_col <- annotation_col[desired_order, , drop = FALSE]
annotation_row <- annotation_row[desired_order_row, , drop = FALSE]
if (!identical(rownames(annotation_col), colnames(heatmap_data))) stop("Column annotation order does not match heatmap data")
if (!identical(rownames(annotation_row), rownames(heatmap_data))) stop("Row annotation order does not match heatmap data")

order_palette <- c(
  "Anguilliformes" = "#e0f3db", "Carangiformes" = "#ccebc5",
  "Mugiliformes" = "#fff7ec", "Aulopiformes" = "#fee8c8",
  "Lutjaniformes" = "#fdd49e", "Pleuronectiformes" = "#fff7fb",
  "Gadiformes" = "#ece7f2", "Clupeiformes" = "#d0d1e6",
  "Spariformes" = "#fde0dd", "Siluriformes" = "#fcc5c0",
  "Chaetodontiformes" = "#f9f9ed", "Perciformes" = "#fcfbde",
  "Carcharhiniformes" = "#fffacc", "Gobiiformes" = "#d9f0a3",
  "Tetraodontiformes" = "#f2e0df", "Eupercaria incertae sedis" = "#e5f5f9",
  "Scombriformes" = "#ccece6", "Myctophiformes" = "#e0ecf4",
  "Carangaria incertae sedis" = "#bfd3e6", "Myliobatiformes" = "#f7fcf0"
)
order_palette <- order_palette[unique(as.character(annotation_col$Order))]
type_palette <- c("Actinopteri" = "#b9a8c2")
family_palette <- c(
  "unclassified Reovirales" = "#f4cccc", "Nodaviridae" = "#e69138",
  "Iridoviridae" = "#f6b26b", "Vibrionaceae" = "#f768a1",
  "Moraxellaceae" = "#b2e3f1", "Pseudomonadaceae" = "#b4a7d6",
  "Glugeidae" = "#fce5cd"
)
family_palette <- family_palette[unique(as.character(annotation_row$Family))]
pathogen_palette <- c("RNA virus" = "#fbb4ae", "DNA virus" = "#ffd9a8",
                      "Bacteria" = "#b3cde4", "Eukaryota" = "#cceac4")

annotation_col$Order <- factor(annotation_col$Order, levels = names(order_palette))
annotation_col$Type <- factor(annotation_col$Type, levels = names(type_palette))
annotation_row$Family <- factor(annotation_row$Family, levels = names(family_palette))
annotation_row$Viruses <- factor(annotation_row$Viruses, levels = names(pathogen_palette))

col_run <- rle(as.character(annotation_col$Order))
gaps_col <- cumsum(col_run$lengths)
if (length(gaps_col) > 1) gaps_col <- gaps_col[-length(gaps_col)] else gaps_col <- NULL
row_run <- rle(as.character(annotation_row$Family))
gaps_row <- cumsum(row_run$lengths)
if (length(gaps_row) > 1) gaps_row <- gaps_row[-length(gaps_row)] else gaps_row <- NULL

plot_obj <- pheatmap(
  heatmap_data,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = list(Order = order_palette, Type = type_palette,
                           Family = family_palette, Viruses = pathogen_palette),
  color = colorRampPalette(c("#f7f7f7", "#b80d57"))(200),
  breaks = seq(0, 2.1, length.out = 201),
  legend_breaks = c(0, 0.5, 1, 1.5, 2),
  legend_labels = c("0", "0.5", "1", "1.5", "2"),
  labels_col = gsub("_", " ", colnames(heatmap_data)),
  labels_row = gsub("_", " ", rownames(heatmap_data)),
  show_colnames = TRUE, show_rownames = TRUE,
  annotation_names_row = FALSE,
  cluster_cols = FALSE, cluster_rows = FALSE,
  gaps_row = gaps_row, gaps_col = gaps_col,
  fontsize = 9, fontsize_row = 8.2, fontsize_col = 7.2,
  angle_col = 90, legend = TRUE, border_color = NA,
  main = NA, silent = TRUE
)

pdf_file <- file.path(script_dir, "Figure6C_heatmap.pdf")
grDevices::pdf(pdf_file, width = 10, height = 10, useDingbats = FALSE)
grid.newpage(); grid.draw(plot_obj$gtable); grDevices::dev.off()
png_file <- file.path(script_dir, "Figure6C_heatmap.png")
grDevices::png(png_file, width = 10 * 600, height = 10 * 600, res = 600,
               units = "px", type = "cairo", bg = "white")
grid.newpage(); grid.draw(plot_obj$gtable); grDevices::dev.off()

write.table(data.frame(metric = c("host_columns", "pathogen_rows", "column_group_gaps",
                                  "row_group_gaps", "first_host_column_retained"),
                       value = c(ncol(heatmap_data), nrow(heatmap_data),
                                 paste(gaps_col, collapse = ";"), paste(gaps_row, collapse = ";"), "TRUE")),
            file.path(script_dir, "Figure6C_layout_metrics.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
