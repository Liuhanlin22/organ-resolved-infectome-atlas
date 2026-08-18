#!/usr/bin/env Rscript

# Figure 4A — organ composition of each pathogen group
# Input:  figure4A.xlsx
# Output: Figure4A.pdf, Figure4A.png, Figure4A_group_percentages.tsv

required_packages <- c("openxlsx", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    "\nInstall them with: install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
})

# ---------------------------- adjustable parameters -------------------------
font_family <- "Arial"

# Legend order and default within-column order.
organ_order <- c(
  "Brain", "Gill", "Intestine", "Liver",
  "Skin", "Muscles", "Spleen", "Kidney"
)

organ_colors <- c(
  "Brain"     = "#07689F",
  "Gill"      = "#88304E",
  "Intestine" = "#74B49B",
  "Liver"     = "#8A79AF",
  "Skin"      = "#D4A5A5",
  "Muscles"   = "#FDA77F",
  "Spleen"    = "#D6604D",
  "Kidney"    = "#649DAD"
)

# Group blocks are arranged by the organ with the largest percentage.
# Within every block, that dominant-organ percentage decreases left to right.
dominant_organ_priority <- c(
  "Brain", "Gill", "Intestine", "Skin",
  "Liver", "Muscles", "Spleen", "Kidney"
)

# Anchor the dominant organ to a common edge within its group block.
# This makes the left-to-right decrease directly comparable.
bottom_anchored_organs <- c("Brain", "Gill", "Intestine", "Liver")
top_anchored_organs <- c("Skin", "Muscles", "Spleen", "Kidney")

figure_width_in <- 14
figure_height_in <- 3
group_label_angle <- 42
group_label_size <- 8.5
axis_title_size <- 10
legend_text_size <- 9
bar_width <- 0.84

# ---------------------------- path handling ---------------------------------
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
setwd(script_dir)
input_file <- file.path(script_dir, "figure4A.xlsx")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

normalize_header <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", " ", x)
  trimws(gsub("[[:space:]]+", " ", x))
}

# ---------------------------- read and validate ------------------------------
dat <- read.xlsx(input_file, sheet = 1, check.names = FALSE)
original_names <- names(dat)
header_keys <- normalize_header(original_names)

find_column <- function(key) {
  hit <- which(header_keys == normalize_header(key))
  if (length(hit) != 1) {
    stop(
      "Column not found or duplicated: ", key,
      "\nDetected columns: ", paste(original_names, collapse = " | ")
    )
  }
  original_names[hit]
}

group_col <- find_column("Group")
organ_cols <- setNames(vapply(organ_order, find_column, character(1)), organ_order)

dat$Group_for_plot <- trimws(as.character(dat[[group_col]]))
if (anyNA(dat$Group_for_plot) || any(!nzchar(dat$Group_for_plot))) {
  stop("The Group column contains empty values.")
}

for (organ in organ_order) {
  values <- suppressWarnings(as.numeric(dat[[organ_cols[[organ]]]]))
  if (any(!is.finite(values))) {
    stop("Non-numeric or missing values detected in organ column: ", organ)
  }
  if (any(values < 0)) stop("Negative values detected in organ column: ", organ)
  dat[[organ]] <- values
}

# ---------------------------- aggregate and calculate -----------------------
# Sum all pathogens belonging to the same Group.
group_sum <- aggregate(
  dat[, organ_order, drop = FALSE],
  by = list(Group = dat$Group_for_plot),
  FUN = sum,
  na.rm = TRUE
)

group_sum$Total <- rowSums(group_sum[, organ_order, drop = FALSE])
if (any(group_sum$Total <= 0)) {
  stop(
    "These groups have zero total abundance and cannot be converted to percentages: ",
    paste(group_sum$Group[group_sum$Total <= 0], collapse = ", ")
  )
}

group_pct <- group_sum
group_pct[, organ_order] <- 100 * sweep(
  as.matrix(group_sum[, organ_order, drop = FALSE]),
  1,
  group_sum$Total,
  "/"
)

# Resolve ties according to dominant_organ_priority.
group_pct$Dominant_organ <- vapply(seq_len(nrow(group_pct)), function(i) {
  values <- as.numeric(unlist(
    group_pct[i, dominant_organ_priority, drop = FALSE],
    use.names = FALSE
  ))
  dominant_organ_priority[which.max(values)]
}, character(1))

group_pct$Dominant_percentage <- vapply(seq_len(nrow(group_pct)), function(i) {
  group_pct[i, group_pct$Dominant_organ[i]]
}, numeric(1))

group_pct$Dominant_rank <- match(
  group_pct$Dominant_organ,
  dominant_organ_priority
)

# Organ blocks follow the requested priority; percentages decline left to right
# inside each block. Alphabetical Group order is only the final tie-breaker.
sort_index <- order(
  group_pct$Dominant_rank,
  -group_pct$Dominant_percentage,
  group_pct$Group
)
group_pct <- group_pct[sort_index, , drop = FALSE]
group_levels <- group_pct$Group

# Build explicit rectangle coordinates because the stack order is allowed to
# differ between group blocks. The dominant organ is placed at the requested
# common baseline (bottom or top) for every column in that block.
plot_data <- do.call(rbind, lapply(seq_len(nrow(group_pct)), function(i) {
  dominant <- group_pct$Dominant_organ[i]

  if (dominant %in% bottom_anchored_organs) {
    stack_order <- c(dominant, setdiff(organ_order, dominant))
  } else if (dominant %in% top_anchored_organs) {
    stack_order <- c(setdiff(organ_order, dominant), dominant)
  } else {
    stack_order <- organ_order
  }

  values <- as.numeric(unlist(
    group_pct[i, stack_order, drop = FALSE],
    use.names = FALSE
  ))
  cumulative <- cumsum(values)

  data.frame(
    Group = group_pct$Group[i],
    Organ = stack_order,
    Percentage = values,
    Stack_rank = seq_along(stack_order),
    ymin = c(0, head(cumulative, -1)),
    ymax = cumulative,
    stringsAsFactors = FALSE
  )
}))
plot_data$Group <- factor(plot_data$Group, levels = group_levels)
plot_data$Organ <- factor(plot_data$Organ, levels = organ_order)
plot_data$Group_index <- match(as.character(plot_data$Group), group_levels)

# Save the exact values and plotting order for verification.
percentage_output <- group_pct[, c(
  "Group", organ_order, "Total", "Dominant_organ", "Dominant_percentage"
)]
write.table(
  percentage_output,
  file = file.path(script_dir, "Figure4A_group_percentages.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8"
)

# ---------------------------- plot -------------------------------------------
p <- ggplot(plot_data, aes(fill = Organ)) +
  geom_rect(
    aes(
      xmin = Group_index - bar_width / 2,
      xmax = Group_index + bar_width / 2,
      ymin = ymin,
      ymax = ymax
    ),
    color = "white",
    linewidth = 0.18
  ) +
  scale_x_continuous(
    breaks = seq_along(group_levels),
    labels = group_levels,
    expand = expansion(add = c(0.5, 0.5))
  ) +
  scale_fill_manual(
    values = organ_colors,
    breaks = organ_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  # Use a coordinate window instead of scale limits. Percentage sums can be
  # 100 + a tiny floating-point error; scale limits would delete the top tile.
  coord_cartesian(ylim = c(0, 100), expand = FALSE, clip = "on") +
  labs(
    x = NULL,
    y = "Relative abundance (%)",
    fill = "Organs"
  ) +
  theme_classic(base_family = font_family, base_size = 9) +
  theme(
    panel.border = element_rect(color = "#666666", fill = NA, linewidth = 0.45),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.text.x = element_text(
      family = font_family,
      size = group_label_size,
      angle = group_label_angle,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      family = font_family,
      size = 8.5,
      color = "black"
    ),
    axis.title.y = element_text(
      family = font_family,
      size = axis_title_size,
      color = "black",
      margin = margin(r = 5)
    ),
    legend.position = "left",
    legend.direction = "vertical",
    legend.title = element_text(
      family = font_family,
      face = "plain",
      size = 13
    ),
    legend.text = element_text(
      family = font_family,
      size = legend_text_size,
      color = "black"
    ),
    legend.key.height = grid::unit(0.35, "cm"),
    legend.key.width = grid::unit(0.45, "cm"),
    legend.spacing.y = grid::unit(0.04, "cm"),
    plot.margin = margin(t = 6, r = 8, b = 5, l = 4)
  ) +
  guides(fill = guide_legend(title.position = "top", ncol = 1, byrow = TRUE))

ggsave(
  file.path(script_dir, "Figure4A.pdf"),
  p,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(script_dir, "Figure4A.png"),
  p,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  dpi = 600,
  bg = "white"
)

