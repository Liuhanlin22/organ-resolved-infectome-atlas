#!/usr/bin/env Rscript

# Figure 5B — total pathogen abundance among fish orders
# Uses 112 fish with one library from each of eight organs.

required <- c("openxlsx", "ggplot2", "writexl")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing R package(s): ", paste(missing, collapse = ", "),
       "\nInstall with: install.packages(c(",
       paste(sprintf("'%s'", missing), collapse = ", "), "))")
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}

input_file <- file.path(script_dir, "Figure5B_input.xlsx")
if (!file.exists(input_file)) stop("Input not found: ", input_file)

dat <- read.xlsx(input_file, sheet = "fish_summary", check.names = FALSE)
numeric_cols <- c("library_count", "order_n_complete_fish",
                  "total_RPM_8_organs", "log10_total_RPM_plus1")
dat[numeric_cols] <- lapply(dat[numeric_cols], as.numeric)

if (nrow(dat) != 112L) stop("Expected 112 complete fish; found ", nrow(dat))
if (any(dat$library_count != 8L)) stop("Every fish must contain exactly eight organ libraries.")
if (any(!is.finite(dat$log10_total_RPM_plus1))) stop("Non-finite transformed abundance detected.")

dat$order_n_complete_fish <- ave(seq_len(nrow(dat)), dat$order, FUN = length)
plot_data <- dat[dat$order_n_complete_fish >= 6, , drop = FALSE]
plot_data$log10_value <- log10(plot_data$total_RPM_8_organs + 1)

split_data <- split(plot_data, plot_data$order)
group_summary <- do.call(rbind, lapply(names(split_data), function(g) {
  x <- split_data[[g]]$log10_value
  raw <- split_data[[g]]$total_RPM_8_organs
  data.frame(
    order = g,
    n = length(x),
    zero_abundance_n = sum(raw == 0),
    median = median(x),
    Q1 = unname(quantile(x, 0.25)),
    Q3 = unname(quantile(x, 0.75)),
    mean = mean(x),
    stringsAsFactors = FALSE
  )
}))
group_summary <- group_summary[order(-group_summary$median, group_summary$order), ]
order_levels <- group_summary$order
plot_data$order <- factor(plot_data$order, levels = order_levels)

kw <- kruskal.test(log10_value ~ order, data = plot_data)
k <- length(order_levels)
N <- nrow(plot_data)
epsilon_squared <- max(0, unname((kw$statistic - k + 1) / (N - k)))

dunn_test_bh <- function(data) {
  data$rank_value <- rank(data$log10_value, ties.method = "average")
  n_total <- nrow(data)
  ties <- table(data$log10_value)
  tie_correction <- 1 - sum(ties^3 - ties) / (n_total^3 - n_total)
  variance_base <- n_total * (n_total + 1) / 12 * tie_correction
  n_group <- table(data$order)
  mean_rank <- tapply(data$rank_value, data$order, mean)
  pairs <- combn(levels(data$order), 2, simplify = FALSE)

  out <- do.call(rbind, lapply(pairs, function(pair) {
    g1 <- pair[1]
    g2 <- pair[2]
    z <- (mean_rank[g1] - mean_rank[g2]) /
      sqrt(variance_base * (1 / n_group[g1] + 1 / n_group[g2]))
    data.frame(
      group1 = g1,
      group2 = g2,
      n1 = as.integer(n_group[g1]),
      n2 = as.integer(n_group[g2]),
      mean_rank1 = unname(mean_rank[g1]),
      mean_rank2 = unname(mean_rank[g2]),
      z = unname(z),
      p_unadjusted = 2 * pnorm(-abs(z)),
      stringsAsFactors = FALSE
    )
  }))
  out$p_BH <- p.adjust(out$p_unadjusted, method = "BH")
  out$significant_BH_0.05 <- out$p_BH < 0.05
  out
}

dunn_results <- if (kw$p.value < 0.05) dunn_test_bh(plot_data) else {
  data.frame(
    group1 = character(), group2 = character(), n1 = integer(), n2 = integer(),
    mean_rank1 = numeric(), mean_rank2 = numeric(), z = numeric(),
    p_unadjusted = numeric(), p_BH = numeric(), significant_BH_0.05 = logical()
  )
}
sig <- dunn_results[dunn_results$significant_BH_0.05, , drop = FALSE]

format_p <- function(p, prefix = "p") {
  if (p < 0.001) paste0(prefix, " < 0.001") else paste0(prefix, " = ", formatC(p, digits = 3, format = "f"))
}

base_y <- max(plot_data$log10_value)
data_range <- diff(range(plot_data$log10_value))
step_y <- max(0.22, data_range * 0.085)

if (nrow(sig) > 0) {
  sig <- sig[order(sig$p_BH), , drop = FALSE]
  sig$x1 <- match(sig$group1, order_levels)
  sig$x2 <- match(sig$group2, order_levels)
  sig$y <- base_y + seq_len(nrow(sig)) * step_y
  sig$label <- vapply(sig$p_BH, format_p, character(1), prefix = "BH p")
  verticals <- rbind(
    data.frame(x = sig$x1, xend = sig$x1, y = sig$y - step_y * 0.12, yend = sig$y),
    data.frame(x = sig$x2, xend = sig$x2, y = sig$y - step_y * 0.12, yend = sig$y)
  )
} else {
  verticals <- data.frame(x = numeric(), xend = numeric(), y = numeric(), yend = numeric())
}

overall_y <- base_y + (nrow(sig) + 1.4) * step_y
y_upper <- overall_y + step_y * 0.6

old_style_colors <- c(
  "Aulopiformes" = "#c06c84",
  "Clupeiformes" = "#9ee6cf",
  "Eupercaria incertae sedis" = "#7ec0e4",
  "Scombriformes" = "#c9b6e4",
  "Tetraodontiformes" = "#8aae92"
)
missing_colors <- setdiff(order_levels, names(old_style_colors))
if (length(missing_colors)) old_style_colors[missing_colors] <- "#BDBDBD"

p <- ggplot(plot_data, aes(x = order, y = log10_value, fill = order)) +
  geom_boxplot(width = 0.60, alpha = 0.70, outlier.shape = NA,
               linewidth = 0.55, color = "black") +
  geom_jitter(width = 0.16, height = 0, size = 2.8,
              alpha = 0.58, color = "black") +
  scale_fill_manual(values = old_style_colors, drop = FALSE) +
  annotate("text", x = 1, y = overall_y,
           label = paste0("Kruskal–Wallis, ", format_p(kw$p.value)),
           hjust = 0, vjust = 0, family = "Arial", size = 3.5) +
  labs(
    x = "Fish order",
    y = expression(log[10] * "(total pathogen RPM + 1)")
  ) +
  scale_x_discrete(labels = function(x) sub(" incertae sedis$", "\nincertae sedis", x)) +
  scale_y_continuous(limits = c(0, y_upper), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1,
                               color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 11),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    legend.position = "none",
    plot.margin = margin(10, 12, 8, 10)
  )

if (nrow(sig) > 0) {
  p <- p +
    geom_segment(data = sig,
                 aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.45) +
    geom_segment(data = verticals,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, linewidth = 0.45) +
    geom_text(data = sig,
              aes(x = (x1 + x2) / 2, y = y + step_y * 0.05, label = label),
              inherit.aes = FALSE, vjust = 0, family = "Arial", size = 3.0)
}

kw_table <- data.frame(
  test = "Kruskal-Wallis rank sum test",
  statistic_H = unname(kw$statistic),
  df = unname(kw$parameter),
  p_value = kw$p.value,
  epsilon_squared = epsilon_squared,
  N = N,
  number_of_orders = k,
  stringsAsFactors = FALSE
)

stats_file <- file.path(script_dir, "Figure5B_statistics.xlsx")
writexl::write_xlsx(
  list(
    group_summary = group_summary,
    Kruskal_Wallis = kw_table,
    Dunn_BH_all = dunn_results,
    Dunn_BH_significant = sig[, intersect(names(dunn_results), names(sig)), drop = FALSE]
  ),
  stats_file
)

ggsave(file.path(script_dir, "Figure5B.pdf"), p,
       width = 5, height = 6.1, units = "in", device = grDevices::cairo_pdf,
       bg = "white")
ggsave(file.path(script_dir, "Figure5B.png"), p,
       width = 5, height = 6.1, units = "in", dpi = 600, bg = "white")

