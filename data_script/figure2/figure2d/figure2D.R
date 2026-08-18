#!/usr/bin/env Rscript

# Figure 2D - association between Toxocara sp. CTSZ-2 and
# Nanhai rhabdovirus 4 abundance

library(readxl)
library(ggplot2)

# -----------------------------------------------------------------------------
# Input and output
# -----------------------------------------------------------------------------
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
  }

  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    script_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(script_path)) {
      return(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

work_dir <- get_script_dir()
input_file <- file.path(work_dir, "figure2d.xlsx")
output_pdf <- file.path(work_dir, "figure2D.pdf")
output_png <- file.path(work_dir, "figure2D.png")

x_col <- "Toxocara_sp._CTSZ-2"
y_col <- "Nanhai_rhabdovirus_4"

# -----------------------------------------------------------------------------
# Read and validate data
# -----------------------------------------------------------------------------
raw_data <- read_excel(input_file, sheet = 1)

required_cols <- c("lib_id", x_col, y_col)
missing_cols <- setdiff(required_cols, colnames(raw_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

data <- data.frame(
  lib_id = as.character(raw_data[["lib_id"]]),
  x = suppressWarnings(as.numeric(raw_data[[x_col]])),
  y = suppressWarnings(as.numeric(raw_data[[y_col]]))
)
data <- data[complete.cases(data[, c("x", "y")]), , drop = FALSE]

if (nrow(data) < 3) {
  stop("At least three complete observations are required for regression.")
}

# -----------------------------------------------------------------------------
# Linear regression and statistics
# -----------------------------------------------------------------------------
lm_model <- lm(y ~ x, data = data)
model_summary <- summary(lm_model)
r_squared <- unname(model_summary$r.squared)
p_value <- unname(coef(model_summary)["x", "Pr(>|t|)"])

format_p <- function(p) {
  if (p < 0.001) {
    out <- formatC(p, format = "e", digits = 3)
    out <- sub("e-0", "e-", out, fixed = TRUE)
    out <- sub("e+0", "e+", out, fixed = TRUE)
    return(out)
  }
  formatC(p, format = "f", digits = 3)
}

stats_text <- sprintf(
  "R² = %.3f\nP-value = %s",
  r_squared,
  format_p(p_value)
)

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
figure2D <- ggplot(data, aes(x = x, y = y)) +
  geom_smooth(
    method = "lm", formula = y ~ x,
    se = TRUE, level = 0.95,
    color = "#B22222", fill = "#AFAFAF", alpha = 0.50,
    linewidth = 0.9
  ) +
  geom_point(
    shape = 16, size = 3.2,
    color = "#65A9C5", alpha = 0.80
  ) +
  annotate(
    "text",
    x = 0.585, y = 3.10,
    label = stats_text,
    hjust = 0, vjust = 1,
    family = "Arial", size = 5.0,
    lineheight = 1.15
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = "d",
    hjust = 1.55, vjust = 0.25,
    family = "Arial", fontface = "bold", size = 12
  ) +
  scale_x_continuous(
    limits = c(0.50, 2.80),
    breaks = seq(1.0, 2.5, by = 0.5),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(-0.20, 3.30),
    breaks = 0:3,
    expand = expansion(mult = 0)
  ) +
  labs(
    x = expression(paste(
      "Abundance of ", italic("Toxocara sp. CTSZ-2"), ": log(RPM+1)"
    )),
    y = expression(atop(
      paste("Abundance of ", italic("Nanhai rhabdovirus 4")),
      "log(RPM+1)"
    ))
  ) +
  coord_cartesian(clip = "off") +
  theme_grey(base_size = 12, base_family = "Arial") +
  theme(
    plot.title = element_blank(),
    panel.background = element_rect(fill = "#E5E5E5", colour = "black", linewidth = 0.55),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.55),
    panel.grid.major = element_line(colour = "white", linewidth = 0.65),
    panel.grid.minor = element_line(colour = "white", linewidth = 0.40),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 0.45),
    axis.ticks.length = grid::unit(2.2, "pt"),
    axis.text = element_text(colour = "black", size = 11),
    axis.title.x = element_text(colour = "black", size = 12, margin = margin(t = 5)),
    axis.title.y = element_text(colour = "black", size = 12, margin = margin(r = 7)),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(t = 17, r = 10, b = 8, l = 31)
  )

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
ggsave(
  output_pdf, figure2D,
  width = 6.20, height = 4.45, units = "in",
  device = cairo_pdf, bg = "white"
)

ggsave(
  output_png, figure2D,
  width = 6.20, height = 4.45, units = "in",
  dpi = 600, bg = "white"
)
